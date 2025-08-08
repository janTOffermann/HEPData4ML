import os
import glob, uuid, itertools
import numpy as np
import h5py as h5
import ROOT as rt
import uproot as ur
import subprocess as sub
from util.calcs import embed_array
from util.buffer import IndexableLazyLoader
from util.qol_utils.progress_bar import printProgressBarColor
from util.hepmc.hepmc import ExtractHepMCEvents, ExtractHepMCParticles, ParticleToVector, ParticleToProductionVertex, ParticleToEndVertex, IsStable, GetParticleID
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
        self.nevents_per_chunk = 100 # Can be configured. Affects memory footprint.
        self.nparticles_max = int(1e4) # TODO: This is some hardcoded max number of particles to be read in from HepMC. Should be plenty.
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

                self.WriteToDataBuffer(j, '{}.Production.Xmu'.format(self.stable_truth_particle_name), np.vstack([
                    [getattr(vec, method)() for vec in [ParticleToProductionVertex(x) for x in stable_particles]]
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

                    self.WriteToDataBuffer(j, '{}.Production.Xmu'.format(key), np.vstack([
                        [getattr(vec, method)() for vec in [ParticleToProductionVertex(x) for x in truth_selected_particles]]
                        for method in ['T','X','Y','Z']
                    ]).T,
                                        dimensions={1:self.nparticles_truth_selected}
                    )

                    self.WriteToDataBuffer(j,'{}.Stable'.format(key),[IsStable(x) for x in truth_selected_particles],
                                            dimensions={1:self.nparticles_truth_selected}, dtype=np.dtype('bool')
                    )

                    self.WriteToDataBuffer(j, '{}.Decay.Xmu'.format(key), np.vstack([
                        [getattr(vec, method)() for vec in [ParticleToEndVertex(x) for x in truth_selected_particles]]
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

                            # another opportunity to add multiplicity, if we haven't already
                            self.WriteToDataBuffer(j,'{}.N'.format(delphes_type),len(delphes_t))


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

        #TODO: Perform some special post-processing if things like JetFinder.Leading() were used,
        #      to strip down an unnecessary dimension.


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
