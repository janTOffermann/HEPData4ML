import glob, itertools
import numpy as np
import h5py as h5
import ROOT as rt
from util.math.embedding import embed_array
from util.buffer.input import UprootTreeLoader
from util.buffer.output import OutputBuffer, RootOutputBuffer
from util.qol_utils.progress_bar import printProgressBarColor
from util.hepmc.hepmc import ExtractHepMCEvents, ExtractHepMCParticles, ParticleToProductionVertex, ParticleToEndVertex, ParticleToMomenta
from typing import Union, Optional, List, Any, TYPE_CHECKING
from util.misc.timing import profile_method, profile_block
import util.hdf5.hdf5 as h5util
import util.root.utils as rootutil

if(TYPE_CHECKING):
    import sys
    from util.config.config import Configurator
    from util.hepmc.setup import HepMCSetup
    setup = HepMCSetup(verbose=False)
    python_dir = setup.GetPythonDirectory()
    if(python_dir not in sys.path):
        sys.path = [setup.GetPythonDirectory()] + sys.path # prepend, to make sure we pick this one up first
    from pyHepMC3 import HepMC3 as hm
    from util.metadata.meta import MetaDataHandler

class Processor:
    """
    Convert HepMC or Delphes/ROOT file, representing an event, into an HDF5 file where each entry/event corresponds
    with a single jet. (Can also skip jet clustering entirely).
    """
    def __init__(self, configurator:'Configurator'):
        self.configurator = configurator
        self.print_prefix = 'Processor'
        self.delphes = False # will be set to True if SetDelphesFiles() is called

        # Set the output format.
        self.output_format = None
        self.output_extension = None
        self._set_output_format()

        self.SetProgressBarPrefix('Converting HepMC3 -> N-tuple:')
        self.suffix = 'Complete'
        self.bl = 50
        self.verbose = False

        self.outdir = ''

        self.stable_truth_particle_name = 'StableTruthParticles' # a "special" name for the stable truth particles collection; this is always present

        self.SetPostProcessing()

        # Various integers for buffer size, number of particles read into memory from HepMC, number saved to file, etc.
        self.buffer_size = configurator.GetReconstructionBufferSize()
        self.nparticles_max = int(1e4) # TODO: This is some hardcoded max number of particles to be read in from HepMC. Should be plenty.
        self.nparticles_stable = self.configurator.GetNPars()['n_stable']
        self.nparticles_truth_selected = self.configurator.GetNPars()['n_truth']
        self.n_delphes = self.configurator.GetNPars()['n_delphes']
        self.buffer = None # will be initialized in Process()
        self.nevents = None

        # truth selector
        self.SetParticleSelection()

        # Metadata handling
        self.metadata_handler = None

        # Temporary vector -- for handling coordinate conversions
        self.tmp_vector = rt.Math.PxPyPzEVector()

        # Compression level for HDF5, 0 (least) to 9 (most)
        self.copts = 0

    def _set_output_format(self):
        self.output_format = self.configurator.GetReconstructionOutputFormat()

        if(self.output_format.lower() == 'hdf5'):
            self.output_extension = 'h5'
        elif(self.output_format.lower() == 'root'):
            self.output_extension = 'root'
        else:
            self._print('Error: Output format {} not recognized.'.format(self.output_format))
            return

    def GetOutputExtension(self):
        return self.output_extension

    def SetH5Compression(self,val:int):
        if(val > 9):
            val = 9
        elif(val < 0):
            val = 0
        self.copts = val

    def SetMetadataHandler(self,handler:'MetaDataHandler'):
        self.metadata_handler = handler

    def _identity_selector(self,event:'hm.GenEvent'):
        """
        A placeholder truth particle selector.
        """
        return {'TruthParticles',event.particles()}

    def SetParticleSelection(self):
        if(self.configurator is None):
            self.truth_selection = self._identity_selector # placeholder
        else:
            self.truth_selection = self.configurator.GetParticleSelection()

    def SetBufferSize(self,val:int):
        self.buffer_size = int(val)

    def SetDelphesFiles(self,val:List[str]):
        if(len(val) > 0):
            self.delphes = True
            self.delphes_files = val

    def SetProgressBarPrefix(self,text:str):
        self.prefix_level1 = text

    def SetPostProcessing(self,post_proc=None):
        if(post_proc is None): post_proc = self.configurator.GetPostProcessing()
        self.post_processing = post_proc

    def SetStatsFile(self,filename:str):
        self.stats_file = filename

    def SetVerbosity(self,flag:bool):
        self.verbose = flag

    def SetOutputDirectory(self,outdir:str):
        self.outdir = outdir

    def _create_buffer(self, output_file):
        if(self.output_format == 'hdf5'):
            if(self.buffer_size > self.nentries):
                self.SetBufferSize(self.nentries)
            self.buffer = OutputBuffer(self.buffer_size)
            self.buffer.SetFilename(output_file)
            self.buffer.SetNEvents(self.nentries)
        else:
            self.buffer = RootOutputBuffer()
            self.buffer.SetFilename(output_file)
            self.buffer.SetTreeName(self.configurator.GetReconstructionTreeName())

    def ProcessFull(self, hepmc_files:Union[List[str],str], output_file:Optional[str]=None, verbosity:int=0):
        output_file = self.Process(hepmc_files,output_file,verbosity)
        self.PostProcess(hepmc_files,output_file)
        return output_file

    @profile_method('Processor.Process')
    def Process(self, hepmc_files:Union[List[str],str], output_file:Optional[str]=None, verbosity:int=0):

        # Parse hepmc_files
        if(type(hepmc_files) == list): hepmc_files = ['{}/{}'.format(self.outdir,x) for x in hepmc_files]
        else:
            hepmc_files = '{}/{}'.format(self.outdir,hepmc_files)
            hepmc_files = glob.glob(hepmc_files,recursive=True)

        if(output_file is None):
            # if(type(hepmc_files) == list): output_file = hepmc_files[0]
            # else: output_file = hepmc_files
            # if(self.delphes): output_file =  output_file.replace('*','').replace('.root','.h5')
            # else: output_file =  output_file.replace('*','').replace('.root','.h5') # TODO: Probably a bug
            output_file = 'events.hepdata4ml.{}'.format(self.output_extension)
        output_file = '{}/{}'.format(self.outdir,output_file)

        if(self.delphes):
            # NOTE: It's important that the Delphes files and truth HepMC files line up!
            #       The way we usually do this, it will be guaranteed.
            #       However, if you borrow functions from here you will have to keep this in mind,
            #       as things will go wrong if final_state_fiels and delphes_files aren't sorted
            #       the same way (or are of different lengths).
            delphes_arr,var_map = self.PrepDelphesArrays()
            self.nentries = len(delphes_arr)

        ## Extract truth particle info from the HepMC files.
        hepmc_events, self.nentries = ExtractHepMCEvents(hepmc_files,get_nevents=True)

        ## Extract particles from the HepMC events
        particles = ExtractHepMCParticles(hepmc_events,self.nparticles_max)

        if(self.truth_selection is None):
            self.truth_selection = {}
        truth_selected_event_particles = {}
        for key,selection in self.truth_selection.items():
            truth_selected_event_particles[key] = ExtractHepMCParticles(hepmc_events,self.nparticles_truth_selected,selection)

        # Set up buffer -- set output file, and total number of entries
        self._create_buffer(output_file)

        progress_bar_chunk = int(self.nentries/20)

        if(verbosity == 1): printProgressBarColor(0,self.nentries, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        for i in range(self.nentries):

            # 0) Write some keys that are the same across all events in this chunk.
            # NOTE: After switching to OutputBuffer usage, have to avoid filling the whole buffer
            #       since this is the 1st key, and that might trigger a flush.
            self.WriteToDataBuffer(i,'SignalFlag',self.configurator.GetSignalFlag())
            self.WriteToDataBuffer(i,'Event.Index',i)

            # 1) Save the stable truth-level particles from the event.
            with profile_block('Processor.Process: Truth Stable'):

                with profile_block('Processor.Process: Truth Stable - 1'):
                    status = np.array([x.status() for x in particles[i]])
                with profile_block('Processor.Process: Truth Stable - 2'):
                    stable_particles = list(itertools.compress(particles[i], status == 1))

                with profile_block('Processor.Process: Truth Stable - 3'):
                # Explicitly fetch/compute 4-momentum components.
                    momenta = np.array([ParticleToMomenta(particle) for particle in stable_particles])
                    stable_particle_momenta = momenta[:,0,:]
                    stable_particle_momenta_cyl = momenta[:,1,:]

                with profile_block('Processor.Process: Truth Stable - 4'):
                    self.WriteToDataBuffer(i,'{}.N'.format(self.stable_truth_particle_name),len(stable_particles))

                    self.WriteToDataBuffer(i, '{}.Pmu'.format(self.stable_truth_particle_name),
                                        stable_particle_momenta,
                                        dimensions={1:self.nparticles_stable}
                    )

                    self.WriteToDataBuffer(i, '{}.Pmu_cyl'.format(self.stable_truth_particle_name),
                                        stable_particle_momenta_cyl,
                                        dimensions={1:self.nparticles_stable}
                    )

                    self.WriteToDataBuffer(i,'{}.PdgId'.format(self.stable_truth_particle_name),[x.pid() for x in stable_particles],
                                        dimensions={1:self.nparticles_stable}, dtype=np.dtype('i4')
                    )

                    self.WriteToDataBuffer(i,'{}.HepMC3Index'.format(self.stable_truth_particle_name),[x.id() for x in stable_particles],
                                        dimensions={1:self.nparticles_stable}, dtype=np.dtype('i4')
                    )

                with profile_block('Processor.Process: Truth Stable - 5'):
                    prod_vertices = np.array([ParticleToProductionVertex(x) for x in stable_particles])

                    self.WriteToDataBuffer(i, '{}.Production.Xmu'.format(self.stable_truth_particle_name),
                                        prod_vertices,
                                        dimensions={1:self.nparticles_stable}
                    )

            # 2) Extract the filtered truth record from the events.
            with profile_block('Processor.Process: Truth Selection'):
                for key,truth_selected_particles in truth_selected_event_particles.items():

                    momenta = np.array([ParticleToMomenta(particle) for particle in truth_selected_particles[i]])
                    truth_particle_momenta = momenta[:,0,:]
                    truth_particle_momenta_cyl = momenta[:,1,:]

                    self.WriteToDataBuffer(i,'{}.N'.format(key),len(truth_selected_particles[i]))

                    self.WriteToDataBuffer(i, '{}.Pmu'.format(key),
                                            truth_particle_momenta,
                                            dimensions={1:self.nparticles_truth_selected}
                    )

                    self.WriteToDataBuffer(i, '{}.Pmu_cyl'.format(key),
                                            truth_particle_momenta_cyl,
                                            dimensions={1:self.nparticles_truth_selected}
                    )

                    self.WriteToDataBuffer(i,'{}.PdgId'.format(key),[x.pid() for x in truth_selected_particles[i]],
                                            dimensions={1:self.nparticles_truth_selected}, dtype=np.dtype('i4')
                    )

                    self.WriteToDataBuffer(i,'{}.HepMC3Index'.format(key),[x.id() for x in truth_selected_particles[i]],
                                            dimensions={1:self.nparticles_truth_selected}, dtype=np.dtype('i4')
                    )

                    prod_vertices = np.array([ParticleToProductionVertex(x) for x in truth_selected_particles[i]])

                    self.WriteToDataBuffer(i, '{}.Production.Xmu'.format(key),
                                        prod_vertices,
                                        dimensions={1:self.nparticles_truth_selected}
                    )

                    self.WriteToDataBuffer(i,'{}.Stable'.format(key),[x.status()==1 for x in truth_selected_particles[i]],
                                            dimensions={1:self.nparticles_truth_selected}, dtype=np.dtype('bool')
                    )

                    end_vertices = np.array([ParticleToEndVertex(x) for x in truth_selected_particles[i]])

                    self.WriteToDataBuffer(i, '{}.Decay.Xmu'.format(key),
                                        end_vertices,
                                        dimensions={1:self.nparticles_truth_selected}
                    )

            # 3) If Delphes was run, we will also extract the relevant information.
            #    Note that PrepDelphesArrays() has been called earlier, if delphes=True.
            #    That's where the Delphes ROOT files are already read and prepared for access
            #    via uproot.
            #
            #    Note that we choose to store these objects as four-momenta, rather than explicitly
            #    storing their components. This is possibly a relevant detail as some objects have
            #    attributes such as pt, while others have Et. These are the same *if* we assume
            #    the objects themselves to be massless.
            if(self.delphes):
                with profile_block('Processor.Process: Delphes'):

                    for k,delphes_type in enumerate(var_map.keys()): # loop over different kinds of Delphes collections
                        is_track = False

                        with profile_block('Processor.Process: Delphes - {}'.format(delphes_type)):

                            if('missinget' in delphes_type.lower()):
                                self.n_delphes[k] = 1 # TODO: Would be nice to eliminate this dimension altogether

                            # Not all objects have all fields, so we do a lot of checking here.
                            if('pt' in var_map[delphes_type].keys()):

                                delphes_pt  = delphes_arr[var_map[delphes_type]['pt' ]][i].to_numpy().astype(float)
                                delphes_eta = delphes_arr[var_map[delphes_type]['eta']][i].to_numpy().astype(float)
                                delphes_phi = delphes_arr[var_map[delphes_type]['phi']][i].to_numpy().astype(float)
                                delphes_m   = np.zeros(delphes_pt.shape)

                                # Rather than use rt.Math.PtEtaPhiMVector, vectorize operations with numpy.
                                # This should be faster (although it's typically nicer to use the ROOT objects to safely
                                # handle the coordinate conversions!). - Jan
                                delphes_px = delphes_pt * np.cos(delphes_phi)
                                delphes_py = delphes_pt * np.sin(delphes_phi)
                                delphes_pz = delphes_pt * np.sinh(delphes_eta)
                                delphes_e  = np.sqrt(np.square(delphes_px) + np.square(delphes_py) + np.square(delphes_pz)) # masses set to zero -> can leave out

                                self.WriteToDataBuffer(i,'{}.N'.format(delphes_type),len(delphes_pt))

                                self.WriteToDataBuffer(i, '{}.Pmu'.format(delphes_type),
                                                    np.column_stack([delphes_e,delphes_px,delphes_py,delphes_pz]),
                                                    dimensions={1:self.n_delphes[k]}
                                )

                                self.WriteToDataBuffer(i, '{}.Pmu_cyl'.format(delphes_type),
                                                    np.column_stack([delphes_pt,delphes_eta,delphes_phi,delphes_m]),
                                                    dimensions={1:self.n_delphes[k]}
                                )

                            if('d0' in var_map[delphes_type].keys()):
                                delphes_d0  = delphes_arr[var_map[delphes_type]['d0']][i].to_numpy()
                                delphes_z0  = delphes_arr[var_map[delphes_type]['z0']][i].to_numpy()
                                delphes_d0e  = delphes_arr[var_map[delphes_type]['errord0']][i].to_numpy()
                                delphes_z0e  = delphes_arr[var_map[delphes_type]['errorz0']][i].to_numpy()

                                self.WriteToDataBuffer(i, '{}.D0'.format(delphes_type), delphes_d0, dimensions={1:self.n_delphes[k]})
                                self.WriteToDataBuffer(i, '{}.D0.Error'.format(delphes_type), delphes_d0e, dimensions={1:self.n_delphes[k]})
                                self.WriteToDataBuffer(i, '{}.Z0'.format(delphes_type), delphes_z0, dimensions={1:self.n_delphes[k]})
                                self.WriteToDataBuffer(i, '{}.Z0.Error'.format(delphes_type), delphes_z0e, dimensions={1:self.n_delphes[k]})

                            if('xd' in var_map[delphes_type].keys()):
                                delphes_xd  = delphes_arr[var_map[delphes_type]['xd']][i].to_numpy()
                                delphes_yd  = delphes_arr[var_map[delphes_type]['yd']][i].to_numpy()
                                delphes_zd  = delphes_arr[var_map[delphes_type]['zd']][i].to_numpy()

                                # store 3-position of closest approach as a vector (Xd, Yd, Zd). Unfortunately Delphes' ParticlePropagator computes Td but doesn't save it...?!
                                #  NOTE: Could consider adding in Td on my own branch of Delphes -- already use this for some other things.
                                self.WriteToDataBuffer(i, '{}.Xdi'.format(delphes_type),
                                                    np.column_stack([delphes_xd, delphes_yd, delphes_zd]),
                                                    dimensions={1:self.n_delphes[k]}
                                )
                                is_track = True # only tracks have this component

                            # In principle, d0, dz and phi give a different way to get Xdi.
                            # TODO: Double-check this!
                            elif('d0' in var_map[delphes_type].keys() and 'z0' in var_map[delphes_type].keys() and 'phi' in var_map[delphes_type].keys()):
                                # d0, z0 and phi already extracted above
                                delphes_xd = delphes_d0 * np.cos(delphes_phi)
                                delphes_yd = delphes_d0 * np.sin(delphes_phi)
                                delphes_zd = delphes_z0
                                self.WriteToDataBuffer(i, '{}.Xdi'.format(delphes_type),
                                                    np.column_stack([delphes_xd, delphes_yd, delphes_zd]),
                                                    dimensions={1:self.n_delphes[k]}
                                )

                            if('charge' in var_map[delphes_type].keys()):
                                delphes_charge = delphes_arr[var_map[delphes_type]['charge']][i].to_numpy()
                                self.WriteToDataBuffer(i, '{}.Charge'.format(delphes_type), delphes_charge, dimensions={1:self.n_delphes[k]}, dtype=float)

                            if('pid' in var_map[delphes_type].keys()):
                                delphes_pid  = delphes_arr[var_map[delphes_type]['pid']][i].to_numpy()
                                self.WriteToDataBuffer(i, '{}.PdgId'.format(delphes_type), delphes_pid, dimensions={1:self.n_delphes[k]}, dtype=np.dtype('i4'))

                            if('eem' in var_map[delphes_type].keys()):
                                delphes_e_em   = delphes_arr[var_map[delphes_type]['eem' ]][i].to_numpy()
                                self.WriteToDataBuffer(i, '{}.E.EM'.format(delphes_type), delphes_e_em, dimensions={1:self.n_delphes[k]}, dtype=float)

                            if('ehad' in var_map[delphes_type].keys()):
                                delphes_e_had  = delphes_arr[var_map[delphes_type]['ehad']][i].to_numpy()
                                self.WriteToDataBuffer(i, '{}.E.Hadronic'.format(delphes_type), delphes_e_had, dimensions={1:self.n_delphes[k]}, dtype=float)

                            if('etrk' in var_map[delphes_type].keys()):
                                delphes_e_trk  = delphes_arr[var_map[delphes_type]['etrk']][i].to_numpy()
                                self.WriteToDataBuffer(i, '{}.E.Track'.format(delphes_type), delphes_e_trk, dimensions={1:self.n_delphes[k]}, dtype=float)

                            # Calorimeter towers indicate their edges in (eta,phi).
                            if('edges' in var_map[delphes_type].keys()):
                                delphes_edges  = delphes_arr[var_map[delphes_type]['edges']][i].to_numpy()
                                # separate eta and phi edges -- I think this is clearer for later reference
                                self.WriteToDataBuffer(i, '{}.Edges.Eta'.format(delphes_type), delphes_edges[:,:2], dimensions={1:self.n_delphes[k]})
                                self.WriteToDataBuffer(i, '{}.Edges.Phi'.format(delphes_type), delphes_edges[:,2:4], dimensions={1:self.n_delphes[k]})

                            # Certain objects record their position in (t,x,y,z). Note that tracks *do not* do this (those are all zero for them).
                            if('x' in var_map[delphes_type].keys() and not is_track):
                                delphes_t  = delphes_arr[var_map[delphes_type]['t' ]][i].to_numpy()
                                delphes_x  = delphes_arr[var_map[delphes_type]['x' ]][i].to_numpy()
                                delphes_y  = delphes_arr[var_map[delphes_type]['y' ]][i].to_numpy()
                                delphes_z  = delphes_arr[var_map[delphes_type]['z' ]][i].to_numpy()

                                self.WriteToDataBuffer(i, '{}.Xmu'.format(delphes_type),
                                                    np.column_stack([delphes_t,delphes_x,delphes_y,delphes_z]),
                                                    dimensions={1:self.n_delphes[k]}
                                )

                                # another opportunity to add multiplicity, if we haven't already
                                self.WriteToDataBuffer(i,'{}.N'.format(delphes_type),len(delphes_t))

            if(verbosity == 1):
                if(((i+1)%progress_bar_chunk == 0) or i+1==self.nentries):
                    printProgressBarColor(i+1,self.nentries, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        # Final flush, in case there are any stragglers in the buffer.
        self.buffer.flush()

        # Close
        self.buffer.close()

        return output_file

    def AddKeyToDataBuffer(self,key:str,value:Union[int,float,np.ndarray,list],dtype:Optional[Union[str,np.dtype]]=None,dimensions:dict=None):
        if(key in self.buffer.keys()):
            return
        value_array = np.asarray(value)

        if(dtype is None):
            dtype = value_array.dtype

        # Create the buffer shape.
        # Note that the 1st dimension, the "event index" dimension,
        # is the buffer size. That will be handled internally by
        # OutputBuffer, so we just need to give the *rest*
        # of the shape.

        if value_array.shape == ():  # Scalar value
            # In buffer, will create shape (N,) where N = buffer size
            buffer_shape = value_array.shape
        else:
            # The user can optionally specify dimensions via a dictionary,
            # otherwise they are inferred.
            # NOTE: Keep in mind that the 1st dimension here is a dummy dimension,
            #       just keeping it to make the code a bit more readable.
            buffer_shape = list((1,) + value_array.shape)
            if(dimensions is not None):
                for idx,val in dimensions.items():
                    try:
                        buffer_shape[idx] = int(val)
                    except:
                        pass # TODO: Add warning
            buffer_shape = tuple(buffer_shape[1:]) # remove dummy dimension, OutputBuffer will prepend buffer size dimension

        self.buffer.create_array(key,buffer_shape,dtype)
        return

    @profile_method('Processor.WriteToDataBuffer')
    def WriteToDataBuffer(self,event_index:Optional[int],key:str,value:Union[int,float,np.ndarray,list],dtype:Optional[Union[str,np.dtype]]=None,dimensions:dict=None):
        if(key not in self.buffer.keys()):
            self.AddKeyToDataBuffer(key,value,dtype,dimensions)

        self.buffer.set(key,event_index,value)
        return

    @profile_method('Processor.PostProcess')
    def PostProcess(self,hepmc_file:str, ntuple_file:str, prepend_output_directory=False):
        if(self.post_processing is None):
            return
        for post_proc in self.post_processing:
            if(post_proc is None): continue
            post_proc.SetConfigurator(self.configurator)
            post_proc.SetMetadataHandler(self.metadata_handler)
            if(prepend_output_directory): # TODO: Clean this up? A little inconsistent
                ntuple_file = '{}/{}'.format(self.outdir,ntuple_file)
                hepmc_file = '{}/{}'.format(self.outdir,hepmc_file)
            post_proc(hepmc_file,ntuple_file,ntuple_file)
        return

    @profile_method('Processor.PrepDelphesArrays')
    def PrepDelphesArrays(self,):
        types = self.configurator.GetDelphesObjects()
        delphes_object_components = [
            'PT','Eta','Phi','ET', 'MET', # momentum componenets
            'D0','ErrorD0','DZ','ErrorDZ', # impact parameters and associated uncertainties (for tracks). NOTE: Delphes seems to have misspelt "Z0" -> "DZ"!
            'X', 'Y', 'Z', 'T', # 4-position -- relevant for non-track objects (tracks parameterized differently)
            'Xd', 'Yd', 'Zd', # 3-position of track position of closest approach to z-axis
            'Eem','Ehad','Etrk', # energy depositions -- relevant for EFlow objects (possibly quite detector card-specific!)
            'Edges[4]', # edges in eta and phi, for calorimeter towers -- format is (etaMin, etaMax, phiMin, phiMax)
            'Charge', 'PID'
        ]
        delphes_keys = ['{x}.{y}'.format(x=x,y=y) for x in types for y in delphes_object_components]
        delphes_tree = 'Delphes'
        delphes_files = ['{}/{}'.format(self.outdir, x) for x in self.delphes_files]

        delphes_arr = UprootTreeLoader(delphes_files, delphes_tree, delphes_keys) # TODO: Switch to using some kind of buffering. This class loads everyhing.
        delphes_keys = delphes_arr.fields # keeps only the fields that actually exist!

        # Create var_map as before
        var_map = {key:{} for key in types}
        for branch in delphes_keys:
            if '.' in branch:
                key, var = branch.split('.')
                var = var.lower()
                if(var=='et' or var=='met'):
                    var = 'pt'
                elif(var=='dz'):
                    var = 'z0'
                elif(var=='errordz'):
                    var = 'errorz0'
                elif(var=='edges[4]'):
                    var = 'edges'
                var_map[key][var.lower()] = branch
        return delphes_arr, var_map

    ################################################################
    # Some member functions that will behave like utilities,
    # i.e. we'll call them separate from the main Processor routine.
    ################################################################

    def ConcatenateNtuples(self,input_files:list,output_file):
            """
            This is a utility function for concatenating output ntuples.
            """
            output_filepath = '/'.join((self.outdir,output_file))
            if(self.output_format.lower() == 'hdf5'):
                compression_opts = 1 # TODO fetch from configuration?
                h5util.ConcatenateH5(input_files,output_filepath,copts=compression_opts,delete_inputs=True,ignore_keys=['Event.Index'],verbose=False,silent_drop=True)
            elif(self.output_format.lower() == 'root'):
                rootutil.ConcatenateRootTreeFiles(input_files,output_filepath,self.buffer.GetTreeName(),['Event.Index'])
            else:
                self._print('Error: ConcatenateNtuples() not implemented for format {}.'.format(self.output_format))
                return

    def AddEventIndices(self,input_file,key='Event.Index',offset=0):
        key = 'Event.Index' # TODO: Make this member variable?
        if(self.output_format.lower() == 'hdf5'):
            compression_opts = 1 # TODO fetch from configuration?
            h5util.AddEventIndices(input_file,cwd=self.outdir,copts=compression_opts,key=key,offset=offset)
        elif(self.output_format.lower() == 'root'):
            rootutil.AddEventIndices(input_file,self.buffer.GetTreeName(),cwd=self.outdir,index_branch_name='Event.Index',offset=offset)
        else:
            self._print('Error: AddEventIndices() not implemented for format {}.'.format(self.output_format))
            return

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return
