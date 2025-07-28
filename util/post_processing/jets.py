# The purpose of this code is to apply the Johns Hopkins top tagger (arXiv:0806.0848 [hep-ph])
# to the jets in the dataset.
import json
import numpy as np
import h5py as h5
from typing import Any, Optional, List # experimenting with typing
from util.fastjet.jetfinderbase import JetFinderBase
from util.qol_utils.progress_bar import printProgressBarColor
from util.buffer import Buffer
from util.meta import AddMetaDataWithReference, GetMetadata

import util.post_processing.utils.ghost_association as ghost_assoc
import util.post_processing.utils.softdrop as softdrop
import util.post_processing.utils.jhtagger as jhtagger
import util.post_processing.utils.jet_filter as jet_filter
import util.post_processing.utils.containment as containment
import util.post_processing.utils.simple_btag as simple_btag

class JetFinder(JetFinderBase):
    """
    This class uses the Fastjet library to perform jet clustering, on some
    (arbitrary) set of inputs representing four-momenta of some objects.
    """

    def __init__(self, input_collections:List[str]=['StableTruthParticles'], jet_algorithm:str='anti_kt',radius:float=0.4, jet_name:str='AK04Jets', n_jets_max:int=10,save_constituents:bool=True, fastjet_dir:Optional[str]=None,verbose:bool=False):

        self.status = False

        self.SetInputCollections(input_collections)
        self.jet_algorithm_name = jet_algorithm
        self.jet_name = jet_name
        self.radius = radius
        self.constituents_flag = save_constituents
        self.n_jets_max = n_jets_max # max number of jets to save per event (will be pt-ordered)
        self.n_constituents_max = 200 # max number of constituents to save per jet

        self.fastjet_dir = fastjet_dir
        self.fastjet_init_flag = False

        self.buffer = Buffer(100) # TODO: Make buffer size configurable

        self.input_collection_arrays = None
        self.input_vecs = None
        self.constituent_indices_dict = None

        self.SetVerbosity(verbose)
        self.configurator = None

        self.print_prefix = '\n\tJetFinder'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tRunning JetFinder:'
        self.progress_bar_suffix = 'Complete'
        self.setup = None
        self.tagger = None

        self.copts = 9

        self.error = False

        self._i = 0
        self.processors = [] # supposedly this is an example of an "observer pattern"

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return

    def SetVerbosity(self,flag:bool):
        self.verbose = flag

    def SetH5EventFile(self,file:str):
        self.h5_file = file
        self.buffer.SetFilename(self.h5_file)

    def SetInputCollections(self,collections:List[str]):
        if(type(collections) != list):
            collections = [collections]
        collections = ['{}.Pmu'.format(collection) for collection in collections]

        self.input_collections = collections

    def SetNConstituentsMax(self,n:int):
        self.n_constituents_max = n

    def SetRadius(self,radius:float):
        self.radius = radius

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def SetUserInfo(self,val):
        self.user_info = val

    def AddUserInfo(self,idx,val):
        if(self.user_info is None):
            self.user_info = {}
        self.user_info[idx] = val

    def GetUserInfo(self,idx):
        return self.user_info[idx]

    def ClearUserInfo(self):
        self.user_info = None

    def _input_consistency_check(self):
        f = h5.File(self.h5_file,'r')
        keys = list(f.keys())
        f.close()

        cleaned_collections = []
        for collection in self.input_collections:
            # we will be using the Cartesian versions of each collection for clustering
            found = collection in keys
            if(not found):
                self._print('Warning: Did not find key {} in file {}. Disabling as input...'.format(collection,self.h5_file))
            else:
                cleaned_collections.append(collection)
        self.input_collections = cleaned_collections
        if(len(self.input_collections)==0):
            self._print('Error: No input collections.')
            return False
        return True

    def Initialize(self): # TODO: May want to consider chunking things and using a buffer? Memory usage will scale better for larger files.
        """
        Reads in the input HDF5 file, and places the required arrays in memory.
        """

        # If already initialized, no need to do it again.
        if(self.status):
            return

        # Input consistency check.
        if(not self._input_consistency_check()):
            self.status = False
            return

        # Initialize fastjet.
        self._initialize_fastjet()
        if(not self.fastjet_init_flag):
            self.status = False
            return

        # Initialize the fastjet jet definition
        self._initialize_jet_definition()

        # Read in the input 4-momenta from the input file.
        # TODO: Currently reading things into memory, which is not ideal for large files.
        #       Should ultimately move towards batching things, which will require
        #       keeping the input file open the whole time.
        f = h5.File(self.h5_file,'r')
        self.input_collection_arrays = {
            key:f[key][:] for key in self.input_collections
        }
        self.nevents = f[self.input_collections[0]].shape[0]

        # With nevents defined, we can initialize the buffer.
        self._initializeBuffer()

        # Optional modification of initialize. May be harnessed by some special configurations.
        self._modifyInitialization()

        self.status = True
        f.close()

    def Process(self):
        self.Initialize()

        if(self.verbose):
            self._print('Process: Number of events = {}'.format(self.nevents))
            printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        for self._i in range(self.nevents):

            # Clear the user info, which marks jets in this event.
            # May be harnessed by some special configurations.
            self.ClearUserInfo()

            # Gather the different input collections together, into one array of four-momenta.
            self.input_vecs = np.vstack([self.input_collection_arrays[key][self._i] for key in self.input_collections]) # NOTE: Using self.input_collections_array.keys() can be dangerous, due to modifications/additions to keys by things like GhostAssociation(). Those should not touch self.input_collections, for this reason.

            # Optional modification of inputs. May be harnessed by some special configurations.
            self._modifyInputs()

            self._clusterJets() # fills self.jets_dict

            # Optional modification of jets. May be harnessed by some special configurations.
            self._modifyJets()

            # Pt-sort the jets, and truncate to fixed length given by self.n_jets_max
            self._ptSort(truncate=True)

            # optionally extract information on jet constituents
            if(self.constituents_flag):
                self._fetchJetConstituents() # fills self.constituent_vectors, self.constituent_indices

                # Optional modification of constituents. May be harnessed by some special configurations.
                self._modifyConstituents()

            # now write to buffer
            # TODO: Will turn the buffer into a more complex object, that outwardly looks like a dictionary
            #       but takes care of chunking/writing itself.
            self._writeToBuffer()

            # Optional extension of writing to buffer. May be harnessed by some special configurations.
            self._modifyWrite()

            if(self.verbose):
                printProgressBarColor(self._i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)
        return

    def _modifyInitialization(self):
        for processor in self.processors:
            processor.ModifyInitialization(self)
        return

    def _modifyInputs(self):
        for processor in self.processors:
            processor.ModifyInputs(self)
        return

    def _modifyJets(self):
        for processor in self.processors:
            processor.ModifyJets(self)
        return

    def _modifyWrite(self):
        for processor in self.processors:
            processor.ModifyWrite(self)

    def _modifyConstituents(self):
        #NOTE: Considering removing this, might cause weird interplay
        #      between processors. Its better to modify the actual
        #      fastjet jet's constituents within _modifyJets(), so that
        #      the processors' handling of the jets isn't interleaved.
        for processor in self.processors:
            processor.ModifyConstituents(self)
        return

    def Flush(self):
        """
        This function simply finishes the writing of our data buffer
        to the output file, by doing a final flush.
        It also writes some metadata to the output file.
        """
        self.buffer.flush()
        self._writeMetadata()
        return

    def _writeMetadata(self):
        key = 'Metadata.JetCollections.InputCollections'
        metadata = GetMetadata(self.h5_file,key=key)
        if(metadata is None):
            metadata = {}
        else:
            metadata = json.loads(metadata[0])
        metadata[self.jet_name] = [x.replace('.Pmu','') for x in self.input_collections]
        # metadata = json.dumps(metadata) # handled within AddMetaDataWithReference
        AddMetaDataWithReference(self.h5_file,value=metadata,key=key,overwrite=True)

    # Using a generic signature -- should consider making the various post-processors inherit from a single parent class!
    def __call__(self,hepmc_file,h5_file,output_file=None,verbose=None, copts=9, key=None):
        if(verbose is not None): self.SetVerbosity(verbose)
        self.SetH5EventFile(h5_file)
        self.Process()
        self.Flush()
        return self.h5_file

    def _initializeBuffer(self):
        self.buffer.SetNEvents(self.nevents)
        self.buffer.create_array('{}.N'.format(self.jet_name),dtype=np.dtype('i4'))
        self.buffer.create_array('{}.Pmu'.format(self.jet_name),shape=(self.n_jets_max,4),dtype=np.dtype('f8'))
        self.buffer.create_array('{}.Pmu_cyl'.format(self.jet_name),shape=(self.n_jets_max,4),dtype=np.dtype('f8'))
        if(self.constituents_flag):
            self.buffer.create_array('{}.Constituents.N'.format(self.jet_name),shape=(self.n_jets_max,),dtype=np.dtype('i4'))
            self.buffer.create_array('{}.Constituents.Pmu'.format(self.jet_name),shape=(self.n_jets_max,self.n_constituents_max,4),dtype=np.dtype('f8'))
            self.buffer.create_array('{}.Constituents.Pmu_cyl'.format(self.jet_name),shape=(self.n_jets_max,self.n_constituents_max,4),dtype=np.dtype('f8'))

            # Also create buffers corresponding to jet constituents' indices w.r.t. the collections they were pulled from.
            # Note that a jet may have used multiple collections -- so we'll keep track of the index of the collection that a constituent came from,
            # as well as its index *within* that collection.
            self.buffer.create_array('{}.Constituents.Collection'.format(self.jet_name),shape=(self.n_jets_max,self.n_constituents_max),dtype=np.dtype('i4'))
            self.buffer.create_array('{}.Constituents.Collection.Index'.format(self.jet_name),shape=(self.n_jets_max,self.n_constituents_max),dtype=np.dtype('i4'))

        return

    def _computeConstituentIndices(self):
        # Precompute collection boundaries once
        n_per_collection = [len(self.input_collection_arrays[key][self._i]) for key in self.input_collections]
        cumulative_lengths = np.cumsum([0] + n_per_collection)

        self.constituent_indices_dict = {}

        for i, jet in self.jets_dict.items():
            raw_indices = [pj.user_index() for pj in jet.constituents()]

            # Vectorized conversion for all indices at once
            raw_indices_array = np.array(raw_indices)
            collection_indices = np.searchsorted(cumulative_lengths[1:], raw_indices_array, side='right')
            local_indices = raw_indices_array - cumulative_lengths[collection_indices]

            # Combine into pairs
            constituent_indices = np.array(list(zip(collection_indices, local_indices)))

            # Store the results
            self.constituent_indices_dict[i] = constituent_indices
        return

    def _writeToBuffer(self,event_index=None):

        if(event_index is None):
            event_index = self._i

        if(len(self.jets_dict) == 0):
            return #TODO: Check that this is OK?

        # Fill jet information in the buffer.
        self.buffer.set('{}.N'.format(self.jet_name),event_index,len(self.jet_vectors))
        self.buffer.set('{}.Pmu'.format(self.jet_name),event_index,np.vstack([self.jet_vectors[i] for i in self.jet_ordering]))
        self.buffer.set('{}.Pmu_cyl'.format(self.jet_name),event_index,np.vstack([self.jet_vectors_cyl[i] for i in self.jet_ordering]))

        # fill jet constituents
        if(self.constituents_flag):

            self.buffer.set('{}.Constituents.N'.format(self.jet_name),event_index,[len(self.constituent_vectors[i]) for i in self.jet_ordering])

            # Figure out the collections and indices of the constituents
            self._computeConstituentIndices()

            # for key,val in self.constituent_indices_dict.items():
            #     print(key,val)

            # Now we loop, as we're embedding what is really jagged information.
            # TODO: Is there another way? I suspect this slows down things a bit.
            for i,j in enumerate(self.jet_ordering):
                self.buffer.set('{}.Constituents.Pmu'.format(self.jet_name),(event_index,i),self.constituent_vectors[j])
                self.buffer.set('{}.Constituents.Pmu_cyl'.format(self.jet_name),(event_index,i),self.constituent_vectors_cyl[j])
                self.buffer.set('{}.Constituents.Collection'.format(self.jet_name),(event_index,i),self.constituent_indices_dict[j][:,0])
                self.buffer.set('{}.Constituents.Collection.Index'.format(self.jet_name),(event_index,i),self.constituent_indices_dict[j][:,1])
        return

    # NOTE: Will define various functions for performing some modifications to clustering or post-processing of results.
    #       Would be nice to put this somewhere else, but I also want to keep the configuration simple.
    #       These will be member functions that return self, so you can do "constructor().function()" instead of just "constructor()"
    #       and in this way chain together a complex configuration without the constructor having to take a huge number of args.

    def PtFilter(self,pt_min=15.):
        self.processors.append(jet_filter.PtFilter(pt_min))
        return self

    def EtaFilter(self,eta_max=2.):
        self.processors.append(jet_filter.EtaFilter(eta_max))
        return self

    def GhostAssociation(self,truth_key,truth_indices,mode='filter',tag_name=None):
        """
        This function performs ghost association, so that
        we will only return jets are that ghost-associated
        to particles from the branch corresponding with truth_key,
        at the indices specified by truth_indices.

        With mode=='filter', it will filter out non-associated jets.
        With mode=='tag', it will save a ghost association flag to the output.

        Returns self, so this can be chained with the constructor.
        """

        self.processors.append(ghost_assoc.GhostAssociator(truth_key,truth_indices,mode,tag_name))
        return self

    def Softdrop(self,z_cut,beta):
        """
        This function performs the softdrop algorithm.

        Returns self, so this can be chained with the constructor.
        """
        self.processors.append(softdrop.Softdrop(z_cut,beta))
        return self

    def IteratedSoftdrop(self,z_cut,beta,dR_cut,max_depth=10,mode='tag'):
        """
        This function performs the iterated softdrop algorithm.

        Returns self, so this can be chained with the constructor.
        """
        self.processors.append(softdrop.IteratedSoftdrop(z_cut,beta,dR_cut,max_depth,mode))
        return self

    def JohnsHopkinsTagger(self,delta_p=0.1,delta_r=0.19,cos_theta_W_max=0.7,top_mass_range=(150.,200.),W_mass_range=(65.,95.), mode='filter',tag_name=None):
        """
        This function performs top-tagging via
        the Johns Hopkins top tagger.

        With mode=='filter', it will filter out non-associated jets.
        With mode=='tag', it will save a ghost association flag to the output.

        Returns self, so this can be chained with the constructor.
        """

        self.processors.append(jhtagger.JohnsHopkinsTagger(delta_p,delta_r,cos_theta_W_max,top_mass_range,W_mass_range,mode,tag_name))
        return self

    def Containment(self,truth_key,truth_indices,delta_r=None,use_rapidity=True, mode='tag',tag_name=None):
        """
        This function performs "containment tagging", by checking
        the DeltaR between the jet and the particle(s) belonging
        to the collection "truth_key", at "truth_indices".

        With mode=='filter', it will filter out non-contained jets.
        With mode=='tag', it will save a containment flag to the output.

        Returns self, so this can be chained with the constructor.
        """
        if delta_r is None:
            delta_r = self.radius

        self.processors.append(containment.ContainmentTagger(truth_key,truth_indices,delta_r,mode,use_rapidity,tag_name))
        return self


    def TrackCountingBTag(self,track_key,track_pt_min=1., delta_r=0.3, track_ip_max=2., sig_min=6.5, ntracks=3, use_3d=False, mode='tag',tag_name=None):
        """
        This function performs a simple b-tagging algorithm,
        taken from Delphes' TrackCountingBTagging module.
        This counts the number of tracks near the jet that
        meet certain criteria on momentum and displacement.

        With mode=='filter', it will filter out non-tagged jets.
        With mode=='tag', it will save a btag flag to the output.

        Returns self, so this can be chained with the constructor.
        """
        if delta_r is None:
            delta_r = self.radius * 0.75

        self.processors.append(simple_btag.TrackCountingBTagging(mode,track_key,track_pt_min,delta_r,track_ip_max,sig_min,ntracks, use_3d, tag_name))
        return self