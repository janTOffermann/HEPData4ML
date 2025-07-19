import glob, uuid, itertools
import numpy as np
import h5py as h5
import ROOT as rt
import uproot as ur
import subprocess as sub
from util.calcs import embed_array
from util.buffer import IndexableLazyLoader
from util.qol_utils.progress_bar import printProgressBarColor
from util.hepmc.hepmc import ExtractHepMCEvents, ExtractHepMCParticles, ParticleToVector, ParticleToProductionVertex, GetParticleID
from typing import Union, Optional, List, TYPE_CHECKING

if(TYPE_CHECKING):
    import sys
    from util.config import Configurator
    from util.hepmc.setup import HepMCSetup
    setup = HepMCSetup(verbose=False)
    setup.PrepHepMC() # will download/install if necessary
    python_dir = setup.GetPythonDirectory()
    if(python_dir not in sys.path):
        sys.path = [setup.GetPythonDirectory()] + sys.path # prepend, to make sure we pick this one up first
    from pyHepMC3 import HepMC3 as hm

class Processor:
    """
    Convert HepMC or Delphes/ROOT file, representing an event, into an HDF5 file where each entry/event corresponds
    with a single jet. (Can also skip jet clustering entirely).
    """
    def __init__(self, configurator:'Configurator'):
        self.configurator = configurator
        self.delphes = False # will be set to True if SetDelphesFiles() is called

        self.SetProgressBarPrefix('Converting HepMC3 -> HDF5:')
        self.suffix = 'Complete'
        self.bl = 50
        self.verbose = False

        self.outdir = ''

        self.stable_truth_particle_name = 'StableTruthParticles' # a "special" name for the stable truth particles collection; this is always present

        self.cluster_sequence = None
        self.jets = None
        self.jets_filtered = None

        self.SetPostProcessing()

        # Data buffer
        self.data = {}

        # Various integers for buffer size, number of particles read into memory from HepMC, number saved to file, etc.
        self.nevents_per_chunk = int(1e1) # Can be configured. Affects memory footprint.
        self.nparticles_max = int(2e3) # TODO: This is some hardcoded max number of particles to be read in from HepMC. Should be plenty.
        self.nparticles_stable = self.configurator.GetNPars()['n_stable']
        self.nparticles_truth_selected = self.configurator.GetNPars()['n_truth']
        self.n_delphes = self.configurator.GetNPars()['n_delphes']

        # truth selector
        self.SetParticleSelection()

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

    def SetNentriesPerChunk(self,val:int):
        self.nevents_per_chunk = int(val)

    def SetDelphesFiles(self,val:List[str]):
        if(len(val) > 0):
            self.delphes = True
            self.delphes_files = val

    def SetProgressBarPrefix(self,text:str):
        self.prefix_level1 = text

    def SetPostProcessing(self,post_proc=None):
        if(post_proc is None): post_proc = self.configurator.GetPostProcessing()
        self.post_processing = post_proc
        if(self.post_processing is None):
            return

        # # Determine whether or not to generate files containing indices
        # # of particles that were both in the truth and final-state selections.
        # # Useful for things like particle daughter tracing.
        # # TODO: There might be a neater way to handle this.
        # self.SetRecordFinalStateIndices(False)
        # for p in post_proc:
        #     if(p is None): continue
        #     if(p.RequiresIndexing()):
        #         self.SetRecordFinalStateIndices(True)
        #         break

    def SetStatsFile(self,filename:str):
        self.stats_file = filename

    def SetVerbosity(self,flag:bool):
        self.verbose = flag

    def SetOutputDirectory(self,outdir:str):
        self.outdir = outdir

    def Process(self, hepmc_files:Union[List[str],str], h5_file:Optional[str]=None, verbosity:int=0):
        if(type(hepmc_files) == list): hepmc_files = ['{}/{}'.format(self.outdir,x) for x in hepmc_files]
        else: hepmc_files = '{}/{}'.format(self.outdir,hepmc_files)

        if(h5_file is None):
            if(type(hepmc_files) == list): h5_file = hepmc_files[0]
            else: h5_file = hepmc_files
            if(self.delphes): h5_file =  h5_file.replace('*','').replace('.root','.h5')
            else: h5_file =  h5_file.replace('*','').replace('.root','.h5') # TODO: Probably a bug
        else:
            h5_file = '{}/{}'.format(self.outdir,h5_file)

        if(type(hepmc_files) == str): hepmc_files = glob.glob(hepmc_files,recursive=True)

        if(self.delphes):
            # NOTE: It's important that the Delphes files and truth HepMC files line up!
            #       The way we usually do this, it will be guaranteed.
            #       However, if you borrow functions from here you will have to keep this in mind,
            #       as things will go wrong if final_state_fiels and delphes_files aren't sorted
            #       the same way (or are of different lengths).
            delphes_arr,var_map = self.PrepDelphesArrays()
            nentries = len(delphes_arr)

        ## Extract final-state truth particle info from the HepMC files.
        hepmc_events, nentries = ExtractHepMCEvents(hepmc_files,get_nevents=True)

        # Some indexing preparation for writing in chunks.
        start_idxs,stop_idxs,ranges = self.PrepIndexRanges(nentries,self.nevents_per_chunk)
        nchunks = len(ranges)

        dsets = None

        if(verbosity == 1): printProgressBarColor(0,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        for i in range(nchunks):
            # Clear the buffer (for safety).
            for key in self.data.keys(): self.data[key][:] = 0

            # 0) Write some keys that are the same across all events in this chunk.
            self.WriteToDataBuffer(None,'SignalFlag',self.configurator.GetSignalFlag())

            # 1) Save the stable truth-level particles from the event.
            # Extract all the particles from the HepMC events, into memory.
            particles = ExtractHepMCParticles(hepmc_events[start_idxs[i]:stop_idxs[i]],self.nparticles_max)

            for j,event_particles in enumerate(particles): # Loop over events in this chunk

                status = np.array([x.status() for x in event_particles])
                stable_particles = list(itertools.compress(event_particles, status == 1))

                # NOTE: Currently supporting both PyHepMC and HepMC, may want to consider
                #       dropping PyHepMC since in practice we won't be switching back and forth.

                stable_particle_vecs = [ParticleToVector(x) for x in stable_particles]
                stable_particle_prod_vertices = [ParticleToProductionVertex(x) for x in stable_particles]

                self.WriteToDataBuffer(j,'{}.N'.format(self.stable_truth_particle_name),len(stable_particle_vecs))

                self.WriteToDataBuffer(j, '{}.Pmu'.format(self.stable_truth_particle_name), np.vstack([
                    [getattr(vec, method)() for vec in stable_particle_vecs]
                    for method in ['E','Px','Py','Pz']
                ]).T,
                                       dimensions={1:self.nparticles_stable}
                )

                self.WriteToDataBuffer(j, '{}.Pmu_cyl'.format(self.stable_truth_particle_name), np.vstack([
                    [getattr(vec, method)() for vec in stable_particle_vecs]
                    for method in ['Pt','Eta','Phi','M']
                ]).T,
                                       dimensions={1:self.nparticles_stable}
                )

                self.WriteToDataBuffer(j,'{}.PdgId'.format(self.stable_truth_particle_name),[GetParticleID(x) for x in stable_particles],
                                       dimensions={1:self.nparticles_stable}, dtype=np.dtype('i4')
                )

                self.WriteToDataBuffer(j, '{}.Xmu'.format(self.stable_truth_particle_name), np.vstack([
                    [getattr(vec, method)() for vec in stable_particle_prod_vertices]
                    for method in ['T','X','Y','Z']
                ]).T,
                                       dimensions={1:self.nparticles_stable}
                )


            # 2) Extract the filtered truth record from the events.
            for key,selection in self.truth_selection.items():

                truth_selected_event_particles = ExtractHepMCParticles(hepmc_events[start_idxs[i]:stop_idxs[i]],self.nparticles_truth_selected,selection)
                for j,truth_selected_particles in enumerate(truth_selected_event_particles): # Loop over events in this chunk

                    # truth_selected_particles = list(itertools.compress(event_particles, truth_selected_status == 1))
                    truth_selected_particle_vecs = [ParticleToVector(x) for x in truth_selected_particles]
                    truth_selected_particle_prod_vertices = [ParticleToProductionVertex(x) for x in truth_selected_particles]

                    self.WriteToDataBuffer(j,'{}.N'.format(key),len(truth_selected_particle_vecs))

                    self.WriteToDataBuffer(j, '{}.Pmu'.format(key), np.vstack([
                        [getattr(vec, method)() for vec in truth_selected_particle_vecs]
                        for method in ['E','Px','Py','Pz']
                    ]).T,
                                        dimensions={1:self.nparticles_truth_selected}
                    )

                    self.WriteToDataBuffer(j, '{}.Pmu_cyl'.format(key), np.vstack([
                        [getattr(vec, method)() for vec in truth_selected_particle_vecs]
                        for method in ['Pt','Eta','Phi','M']
                    ]).T,
                                            dimensions={1:self.nparticles_truth_selected}
                    )

                    self.WriteToDataBuffer(j,'{}.PdgId'.format(key),[GetParticleID(x) for x in truth_selected_particles],
                                            dimensions={1:self.nparticles_truth_selected}, dtype=np.dtype('i4')
                    )

                    self.WriteToDataBuffer(j, '{}.Xmu'.format(key), np.vstack([
                        [getattr(vec, method)() for vec in truth_selected_particle_prod_vertices]
                        for method in ['T','X','Y','Z']
                    ]).T,
                                        dimensions={1:self.nparticles_truth_selected}
                    )

            # 3) Optional: Extract the full event record, store it separately. This is both particles and vertices.

            # 4) If Delphes was run, we will also extract the relevant information.
            #    Note that PrepDelphesArrays() has been called earlier, if delphes=True.
            #    That's where the Delphes ROOT files are already read and prepared for access
            #    via uproot.
            #
            #    Note that we choose to store these objects as four-momenta, rather than explicitly
            #    storing their components. This is possibly a relevant detail as some objects have
            #    attributes such as pt, while others have Et. These are the same *if* we assume
            #    the objects themselves to be massless.
            if(self.delphes):
                for j in range(len(particles)): # TODO: reusing len(particles) (= number of events in chunk), OK but looks kind of hacky
                    for k,delphes_type in enumerate(var_map.keys()): # loop over different kinds of Delphes collections

                        if('missinget' in delphes_type.lower()):
                            self.n_delphes[k] = 1 # TODO: Would be nice to eliminate this dimension altogether

                        # Not all objects have all fields, so we do a lot of checking here.
                        if('pt' in var_map[delphes_type].keys()):

                            delphes_pt  = delphes_arr[var_map[delphes_type]['pt' ]][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_eta = delphes_arr[var_map[delphes_type]['eta']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_phi = delphes_arr[var_map[delphes_type]['phi']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_m   = np.zeros(delphes_pt.shape)

                            delphes_vecs = [rt.Math.PtEtaPhiMVector(*x) for x in zip(delphes_pt,delphes_eta,delphes_phi,delphes_m)]

                            self.WriteToDataBuffer(j,'{}.N'.format(delphes_type),len(delphes_pt))

                            self.WriteToDataBuffer(j, '{}.Pmu'.format(delphes_type), np.vstack([
                                [getattr(vec, method)() for vec in delphes_vecs]
                                for method in ['E','Px','Py','Pz']
                            ]).T,
                                                dimensions={1:self.n_delphes[k]}
                            )

                            self.WriteToDataBuffer(j, '{}.Pmu_cyl'.format(delphes_type), np.vstack([
                                [getattr(vec, method)() for vec in delphes_vecs]
                                for method in ['Pt','Eta','Phi','M']
                            ]).T,
                                                dimensions={1:self.n_delphes[k]}
                            )

                        is_track = False

                        if('d0' in var_map[delphes_type].keys()):
                            delphes_d0  = delphes_arr[var_map[delphes_type]['d0']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_z0  = delphes_arr[var_map[delphes_type]['z0']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_d0e  = delphes_arr[var_map[delphes_type]['errord0']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_z0e  = delphes_arr[var_map[delphes_type]['errorz0']][start_idxs[i]:stop_idxs[i]][j].to_numpy()

                            self.WriteToDataBuffer(j, '{}.D0'.format(delphes_type), delphes_d0, dimensions={1:self.n_delphes[k]})
                            self.WriteToDataBuffer(j, '{}.D0.Error'.format(delphes_type), delphes_d0e, dimensions={1:self.n_delphes[k]})
                            self.WriteToDataBuffer(j, '{}.Z0'.format(delphes_type), delphes_z0, dimensions={1:self.n_delphes[k]})
                            self.WriteToDataBuffer(j, '{}.Z0.Error'.format(delphes_type), delphes_z0e, dimensions={1:self.n_delphes[k]})

                        if('xd' in var_map[delphes_type].keys()):
                            delphes_xd  = delphes_arr[var_map[delphes_type]['xd']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_yd  = delphes_arr[var_map[delphes_type]['yd']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_zd  = delphes_arr[var_map[delphes_type]['zd']][start_idxs[i]:stop_idxs[i]][j].to_numpy()

                            # store 3-position of closest approach as a vector (Xd, Yd, Zd). Unfortunately Delphes' ParticlePropagator computes Td but doesn't save it...?!
                            #  NOTE: Could consider adding in Td on my own branch of Delphes -- already use this for some other things.
                            self.WriteToDataBuffer(j, '{}.Xdi'.format(delphes_type), np.vstack([
                                delphes_xd, delphes_yd, delphes_zd
                            ]).T,
                                                dimensions={1:self.n_delphes[k]}
                            )
                            # self.WriteToDataBuffer(j, '{}.Xd'.format(delphes_type), delphes_xd, dimensions={1:self.n_delphes[k]})
                            # self.WriteToDataBuffer(j, '{}.Yd'.format(delphes_type), delphes_yd, dimensions={1:self.n_delphes[k]})
                            # self.WriteToDataBuffer(j, '{}.Zd'.format(delphes_type), delphes_zd, dimensions={1:self.n_delphes[k]})
                            is_track = True # only tracks have this component

                        # In principle, d0, dz and phi give a different way to get Xdi.
                        # TODO: Double-check this!
                        elif('d0' in var_map[delphes_type].keys() and 'z0' in var_map[delphes_type].keys() and 'phi' in var_map[delphes_type].keys()):
                            # d0, z0 and phi already extracted above
                            delphes_xd = delphes_d0 * np.cos(delphes_phi)
                            delphes_yd = delphes_d0 * np.sin(delphes_phi)
                            delphes_zd = delphes_z0
                            self.WriteToDataBuffer(j, '{}.Xdi'.format(delphes_type), np.vstack([
                                delphes_xd, delphes_yd, delphes_zd
                            ]).T,
                                                dimensions={1:self.n_delphes[k]}
                            )

                        if('charge' in var_map[delphes_type].keys()):
                            delphes_charge = delphes_arr[var_map[delphes_type]['charge']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            self.WriteToDataBuffer(j, '{}.Charge'.format(delphes_type), delphes_charge, dimensions={1:self.n_delphes[k]}, dtype=float)

                        if('pid' in var_map[delphes_type].keys()):
                            delphes_pid  = delphes_arr[var_map[delphes_type]['pid']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            self.WriteToDataBuffer(j, '{}.PdgId'.format(delphes_type), delphes_pid, dimensions={1:self.n_delphes[k]}, dtype=np.dtype('i4'))

                        if('eem' in var_map[delphes_type].keys()): # assume Eem and Ehad together
                            delphes_e_em   = delphes_arr[var_map[delphes_type]['eem' ]][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            self.WriteToDataBuffer(j, '{}.E.EM'.format(delphes_type), delphes_e_em, dimensions={1:self.n_delphes[k]}, dtype=float)

                        if('ehad' in var_map[delphes_type].keys()):
                            delphes_e_had  = delphes_arr[var_map[delphes_type]['ehad']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            self.WriteToDataBuffer(j, '{}.E.Hadronic'.format(delphes_type), delphes_e_had, dimensions={1:self.n_delphes[k]}, dtype=float)

                        if('etrk' in var_map[delphes_type].keys()):
                            delphes_e_trk  = delphes_arr[var_map[delphes_type]['etrk']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            self.WriteToDataBuffer(j, '{}.E.Track'.format(delphes_type), delphes_e_trk, dimensions={1:self.n_delphes[k]}, dtype=float)

                        # Calorimeter towers indicate their edges in (eta,phi).
                        if('edges' in var_map[delphes_type].keys()):
                            delphes_edges  = delphes_arr[var_map[delphes_type]['edges']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            # separate eta and phi edges -- I think this is clearer for later reference
                            self.WriteToDataBuffer(j, '{}.Edges.Eta'.format(delphes_type), delphes_edges[:,:2], dimensions={1:self.n_delphes[k]})
                            self.WriteToDataBuffer(j, '{}.Edges.Phi'.format(delphes_type), delphes_edges[:,2:4], dimensions={1:self.n_delphes[k]})

                        # Certain objects record their position in (t,x,y,z). Note that tracks *do not* do this (those are all zero for them).
                        if('x' in var_map[delphes_type].keys() and not is_track):

                            delphes_t  = delphes_arr[var_map[delphes_type]['t' ]][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_x  = delphes_arr[var_map[delphes_type]['x' ]][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_y  = delphes_arr[var_map[delphes_type]['y' ]][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                            delphes_z  = delphes_arr[var_map[delphes_type]['z' ]][start_idxs[i]:stop_idxs[i]][j].to_numpy()

                            delphes_xvecs = [rt.Math.XYZTVector(*x) for x in zip(delphes_x,delphes_y,delphes_z,delphes_t)]

                            self.WriteToDataBuffer(j, '{}.Xmu'.format(delphes_type), np.vstack([
                                [getattr(vec, method)() for vec in delphes_xvecs]
                                for method in ['T','X','Y','Z']
                            ]).T,
                                                dimensions={1:self.n_delphes[k]}
                            )

            # We have now filled a chunk, time to write it.
            # If this is the first instance of the loop, we will initialize the HDF5 file.
            # NOTE: We assume that after this first loop, we've generated all the necessary keys.
            #       Probably a safe assumption for now.

            if(i == 0):
                dsets = self.PrepH5File(h5_file,nentries,self.data)

            with h5.File(h5_file, 'a') as f:
                for key in dsets.keys():
                    dset = f[key]
                    dset[start_idxs[i]:stop_idxs[i]] = self.data[key][:ranges[i]]
            # if(verbosity == 1): printProgressBarColor(i+1,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)
            if(verbosity == 1): printProgressBarColor(i+1,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        # if(self.diagnostic_plots): self.OutputHistograms()
        return h5_file

    def AddKeyToDataBuffer(self,key:str,value:Union[int,float,np.ndarray,list],dtype:Optional[Union[str,np.dtype]]=None,dimensions:dict=None):
        if(key in self.data.keys()):
            return
        value_array = np.asarray(value)

        if(dtype is None):
            dtype = value_array.dtype

        if value_array.shape == ():  # Scalar value
            # Create shape (N,)
            self.data[key] = np.zeros(self.nevents_per_chunk, dtype=dtype)
        else:
            # Create the buffer shape.
            # The user can optionally specify dimensions via a dictionary,
            # otherwise they are inferred.
            buffer_shape = list((self.nevents_per_chunk,) + value_array.shape)
            if(dimensions is not None):
                for idx,val in dimensions.items():
                    try:
                        buffer_shape[idx] = int(val)
                    except:
                        pass # TODO: Add warning

            # # Collapse buffer_shape to remove any zeros
            # buffer_shape = [x for x in buffer_shape if x!=0]
            self.data[key] = np.zeros(buffer_shape, dtype=dtype)
        return

    def WriteToDataBuffer(self,event_index:Optional[int],key:str,value:Union[int,float,np.ndarray,list],dtype:Optional[Union[str,np.dtype]]=None,dimensions:dict=None):
        if(key not in self.data.keys()):
            self.AddKeyToDataBuffer(key,value,dtype,dimensions)

        value_array = np.asarray(value) # TODO: not sure if needed?
        if(event_index is not None):
            self.data[key][event_index] = embed_array(value_array,self.data[key][event_index].shape)
        else:
            self.data[key][:] = embed_array(value_array,self.data[key].shape)
        return

    def PostProcess(self,hepmc_files:Union[str,List[str]], h5_files:Optional[Union[str,List[str]]]=None):
        if(not isinstance(hepmc_files,list)):
            hepmc_files = [hepmc_files]

        if(self.post_processing is None):
            return
        nfiles = len(hepmc_files)
        # if(truth_files is not None): assert(nfiles == len(truth_files)) # TODO: may have to rework some of this logic
        if(h5_files is not None): assert(nfiles == len(h5_files))

        for post_proc in self.post_processing:
            if(post_proc is None): continue
            # print('\tRunning post processor: {}'.format(post_proc))
            post_proc.SetConfigurator(self.configurator)
            for i in range(nfiles):
                h5_file = None
                if(h5_files is not None): h5_file = '{}/{}'.format(self.outdir,h5_files[i])
                hepmc_file = '{}/{}'.format(self.outdir,hepmc_files[i])
                post_proc(hepmc_file,h5_file,h5_file)
        return

    def MergeEventFilterFlag(self,h5_file,event_filter_flag_files,copts=9):
        if(type(event_filter_flag_files) not in [list,tuple]):
            event_filter_flag_files = [event_filter_flag_files]
        f = h5.File('{}/{}'.format(self.outdir,h5_file),'a')
        keys = list(f.keys())

        try:
            g = [h5.File('{}/{}'.format(self.outdir,x),'r') for x in event_filter_flag_files]
        except:
            pass
            f.close()
            return

        gkeys = list(g[0].keys())

        keys_to_merge = []
        for key in gkeys:
            if(key in keys):
                print('\tWarning: Found key {} among the event_filter_flag keys, but this matches an existing key in the dataset. Skipping.'.format(key))
                continue
            keys_to_merge.append(key)

        for key in keys_to_merge:
            data = np.concatenate([x[key][:] for x in g],axis=0)
            f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)

        f.close()
        for x in g:
            x.close()
        return

    def PrepDelphesArrays(self,):
        types = self.configurator.GetDelphesObjects()
        components = [
            'PT','Eta','Phi','ET', 'MET', # momentum componenets
            'D0','ErrorD0','DZ','ErrorDZ', # impact parameters and associated uncertainties (for tracks). NOTE: Delphes seems to have misspelt "Z0" -> "DZ"!
            'X', 'Y', 'Z', 'T', # 4-position -- relevant for non-track objects (tracks parameterized differently)
            'Xd', 'Yd', 'Zd', # 3-position of track position of closest approach to z-axis
            'Eem','Ehad','Etrk', # energy depositions -- relevant for EFlow objects (possibly quite detector card-specific!)
            'Edges[4]', # edges in eta and phi, for calorimeter towers -- format is (etaMin, etaMax, phiMin, phiMax)
            'Charge', 'PID'
        ]
        delphes_keys = ['{x}.{y}'.format(x=x,y=y) for x in types for y in components]
        delphes_tree = 'Delphes'
        delphes_files = ['{}/{}'.format(self.outdir, x) for x in self.delphes_files]

        delphes_arr = IndexableLazyLoader(delphes_files, delphes_tree, delphes_keys)
        delphes_keys = delphes_arr.fields

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

    def PrepH5File(self,filename,nentries,data_buffer):
        dsets = {}
        with h5.File(filename, 'w') as f:
            for key, val in data_buffer.items():
                shape = list(val.shape)
                shape[0] = nentries
                shape = tuple(shape)
                dsets[key] = f.create_dataset(key, shape, val.dtype,compression='gzip')
        return dsets

    def PrepIndexRanges(self,nentries,nentries_per_chunk):
        nchunks = int(np.ceil(nentries / nentries_per_chunk))
        start_idxs = np.zeros(nchunks,dtype = np.dtype('i8'))
        for i in range(1,start_idxs.shape[0]): start_idxs[i] = start_idxs[i-1] + nentries_per_chunk
        stop_idxs = start_idxs + nentries_per_chunk
        stop_idxs[-1] = nentries
        ranges = stop_idxs - start_idxs
        return start_idxs,stop_idxs,ranges

# Generic function for concatenating HDF5 files.
def ConcatenateH5(input_files,output_file,cwd=None,delete_inputs=False, compression='gzip', copts=9,ignore_keys=None,verbose=False,silent_drop=False):
    if(cwd is not None):
        input_files = ['{}/{}'.format(cwd,x) for x in input_files]
    infiles = [h5.File(f,'r') for f in input_files]
    keys = list(infiles[0].keys())
    if(cwd is not None):
        output_file = '{}/{}'.format(cwd,output_file)
    outfile = h5.File(output_file,'w')

    if(verbose):
        print('\n\t' + 21*"#")
        print("\t### ConcatenateH5 ###")
        for i,infile in enumerate(input_files):
            print('\t  Input {}:\t{}'.format(i,infile))
        print('\tOutput:\t{}'.format(output_file))

    if(ignore_keys is not None):
        keys = [k for k in keys if k not in ignore_keys]
        if(verbose and not silent_drop):
            print('\tExcluding the following keys from concatenation, these will be dropped:')
            for k in ignore_keys:
                print('\t\t{}'.format(k))

    for key in keys:
        data = np.concatenate([f[key] for f in infiles],axis=0)
        outfile.create_dataset(key,data=data,compression=compression,compression_opts=copts)

    outfile.close()

    for f in infiles:
        f.close()

    if(verbose):
        if(delete_inputs):
            print('\tDeleting input files.')
        print('\t' + 21*"#" + '\n')

    if(delete_inputs):
        for f in input_files:
            sub.check_call(['rm',f])

def MergeH5(target_file, input_file, cwd=None, delete_stats_file=False, compression='gzip',copts=9):
    if(cwd is not None):
        target_file = '{}/{}'.format(cwd,target_file)
        input_file = '{}/{}'.format(cwd,input_file)

    f = h5.File(target_file,'a')
    g = h5.File(input_file,'r')
    keys = list(g.keys())
    f_keys = list(f.keys())

    for key in keys:
        if key in f_keys:
            print('Warning: key {} found in {} before merging in {}.'.format(key,target_file,input_file))
            continue
        f.create_dataset(key,data=g[key][:],compression=compression,compression_opts=copts)

    f.close()
    g.close()
    if(delete_stats_file):
        sub.check_call(['rm',input_file])
    return

def RemoveFailedFromHDF5(h5_file, cwd=None, copts=9):
    """
    Remove any failed events from the HDF5 file -- these are identified
    as those with a negative value in the "SignalFlag" dataset.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)

    fname_tmp = h5_file.replace('.h5','_{}.h5'.format(str(uuid.uuid4())))

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    keep_indices = f['SignalFlag'][:] >= 0
    if(np.sum(keep_indices) == N): return

    g = h5.File(fname_tmp,'w')

    for i,key in enumerate(keys):
        g.create_dataset(key, data = f[key][:][keep_indices],compression='gzip',compression_opts=copts)

    f.close()
    g.close()

    sub.check_call(['rm',h5_file])
    sub.check_call(['mv',fname_tmp, h5_file])
    return

# Add a column with event indices.
def AddEventIndices(h5_file,cwd=None,copts=9,key='event_idx',offset=0):
    """
    Add a dataset with event indices to an HDF5 file. This is a sequential index,
    which may be useful if some entries will later be dropped and one wants to
    keep track of which events these were.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    nevents = f['SignalFlag'].shape[0]
    event_indices = np.arange(offset,nevents+offset,dtype=int)
    f.create_dataset(key,data=event_indices, compression='gzip', compression_opts=copts)
    f.close()

def AddConstantValue(h5_file,cwd=None,copts=9,value=0,key='constant_value',dtype=None):
    """
    Add a dataset with some constant value to the HDF5 file. In practice this can
    be used to attach metadata to events -- like the Pythia8 RNG seed that was used.
    This may become useful when multiple datasets are combined.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    nevents = f['SignalFlag'].shape[0]

    if(dtype is None): dtype = np.dtype('f8')
    if(dtype != str):
        data = np.full((nevents),value)
    else:
        data = nevents * [value]
    f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
    f.close()

def AddMetaData(h5_file,cwd=None,value='',key='metadata'):
    """
    Generic function for adding a value to the HDF5 file
    metadata container.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    f.attrs[key] = value
    f.close()

def AddMetaDataWithReference(h5_file,cwd=None,value='',key='metadata',copts=9):
    """
    Adds an entry to the metadata -- if under an existing key, appends it to the list at that key.
    Also creates a column in the dataset that will point to this metadata's index.
    Somewhat redundant for file generation but this type of logic will be useful when concatenating files
    with different entries in the metadata fields.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    nevents = f['SignalFlag'].shape[0]
    metadata = f.attrs
    if(key not in metadata.keys()): f.attrs[key] = [value]
    else: f.attrs[key] = list(f.attrs[key]) + [value] # I think the list <-> array stuff should be OK here
    idx = len(f.attrs[key]) - 1
    f.create_dataset(key,data=np.full(nevents,idx,dtype=np.dtype('i4')),compression='gzip',compression_opts=copts)
    f.close()

def SplitH5(h5_file, split_ratio = (7,2,1), train_name=None, val_name=None, test_name=None, cwd=None, copts=9,verbose=False, seed=0):
    """
    Split an HDF5 file into training, testing and validation samples.
    The split is random (with the RNG seed provided).
    Note that this copies over the full metadata, though with the splitting
    it is technically possible that one of the files has some metadata
    entries with no corresponding events.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    file_dir = '/'.join(h5_file.split('/')[:-1])

    if(train_name is None):
        train_name = 'train.h5'
    if(val_name is None):
        val_name   =   'valid.h5'
    if(test_name is None):
        test_name  =  'test.h5'

    if(cwd is not None):
        train_name = '{}/{}'.format(cwd,train_name)
        val_name = '{}/{}'.format(cwd,val_name)
        test_name = '{}/{}'.format(cwd,test_name)

    if(verbose):
        print("\n\t" + 15*"#")
        print('\t### SplitH5 ###')
        print('\tInput: {}'.format(h5_file))
        print('\tOutputs:')
        for (name,frac) in zip((train_name,val_name,test_name),split_ratio):
            print('\t\t {} \t (fraction of events = {:.2e})'.format(name,frac))
        print('\tSplitting using seed: {}'.format(seed))

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    split_ratio = np.array(split_ratio)
    tot = np.sum(np.array(split_ratio))
    split_ratio = split_ratio / tot
    split_ratio[-1] = 1. - np.sum(split_ratio[:-1])

    names = [train_name, val_name, test_name]
    n = np.array(N * split_ratio,dtype=int)
    diff = np.sum(n) - N
    n[-1] -= diff

    rng = np.random.default_rng(seed)
    index_list = np.array(range(N),dtype=int)
    rng.shuffle(index_list,axis=0)
    indices = []
    start = 0
    stop = 0
    for i in range(len(names)):
        stop += n[i]
        idxs = index_list[start:stop]
        start = stop
        indices.append(idxs)

    for i in range(len(names)):
        g = h5.File(names[i],'w')
        for j,key in enumerate(keys):
                # Putting the data from f into memory (numpy array)
                # since h5py doesn't like non-ordered indices.
                g.create_dataset(key, data=f[key][:][indices[i]],compression='gzip', compression_opts=copts)
        for key in list(f.attrs.keys()):
            g.attrs[key] = f.attrs[key]
        g.close()
    f.close()

    if(verbose):
        print("\t" + 15*"#" + "\n")
    return