# The purpose of this code is to apply the Johns Hopkins top tagger (arXiv:0806.0848 [hep-ph])
# to the jets in the dataset.
import subprocess as sub
import numpy as np
import h5py as h5
# from typing import Dict, Any, Tuple, Optional # experimenting with typing
from util.fastjet import JetFinderBase
from util.qol_utils.progress_bar import printProgressBarColor
from util.calcs import embed_array_inplace
from util.buffer import Buffer

import util.post_processing.utils.ghost_association as ghost_assoc
import util.post_processing.utils.softdrop as softdrop
import util.post_processing.utils.jhtagger as jhtagger
import util.post_processing.utils.jet_filter as jet_filter


class JetFinder(JetFinderBase):
    """
    This class uses the Fastjet library to perform jet clustering, on some
    (arbitrary) set of inputs representing four-momenta of some objects.
    """

    def __init__(self, input_collections=['StableTruthParticles'], jet_algorithm='anti_kt',radius=0.4, jet_name='AK04Jets', n_jets_max=10,save_constituents=True, fastjet_dir=None,verbose=False):

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

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return

    def SetVerbosity(self,flag):
        self.verbose = flag

    def SetH5EventFile(self,file):
        self.h5_file = file
        self.buffer.SetFilename(self.h5_file)

    def SetInputCollections(self,collections):
        if(type(collections) != list):
            collections = [collections]
        collections = ['{}.Pmu'.format(collection) for collection in collections]

        self.input_collections = collections

    def SetNConstituentsMax(self,n):
        self.n_constituents_max = n

    def SetRadius(self,radius):
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
            self.input_vecs = np.vstack([self.input_collection_arrays[key][self._i] for key in self.input_collections]) # NOTE: Using self.input_collections_array.keys() can be dangerous, due to modifications/additions to keys by things like GhostAssociation()

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
        """
        # if(self.verbose):
        self._print('Finishing flushing data to {}.'.format(self.h5_file))
        self.buffer.flush()
        # f = h5.File(self.h5_file,'a')
        # for key,val in self.buffer.items():
        #     d = f.create_dataset(key, data=val, compression='gzip',compression_opts=self.copts)
        # f.close()
        return

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

        # self.buffer['{}.N'.format(self.jet_name)] = np.zeros((self.nevents),dtype=np.dtype('i4'))
        # self.buffer['{}.Pmu'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,4),dtype=np.dtype('f8'))
        # self.buffer['{}.Pmu_cyl'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,4),dtype=np.dtype('f8'))
        # if(self.constituents_flag):
        #     self.buffer['{}.Constituents.N'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max),dtype=np.dtype('i4'))
        #     self.buffer['{}.Constituents.Pmu'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,self.n_constituents_max, 4),dtype=np.dtype('f8'))
        #     self.buffer['{}.Constituents.Pmu_cyl'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,self.n_constituents_max,4),dtype=np.dtype('f8'))
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
            # embed_array_inplace([len(self.constituent_vectors[i]) for i in self.jet_ordering], self.buffer['{}.Constituents.N'.format(self.jet_name)][event_index])

            # Now we loop, as we're embedding what is really jagged information.
            # TODO: Is there another way? I suspect this slows down things a bit.
            for i,j in enumerate(self.jet_ordering):
                self.buffer.set('{}.Constituents.Pmu'.format(self.jet_name),(event_index,i),self.constituent_vectors[j])
                self.buffer.set('{}.Constituents.Pmu_cyl'.format(self.jet_name),(event_index,i),self.constituent_vectors_cyl[j])

                # embed_array_inplace(self.constituent_vectors[j],self.buffer['{}.Constituents.Pmu'.format(self.jet_name)][event_index,i])
                # embed_array_inplace(self.constituent_vectors_cyl[j],self.buffer['{}.Constituents.Pmu_cyl'.format(self.jet_name)][event_index,i])
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