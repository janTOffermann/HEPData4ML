# The purpose of this code is to apply the Johns Hopkins top tagger (arXiv:0806.0848 [hep-ph])
# to the jets in the dataset.
import numpy as np
from typing import Any, Optional, List, Tuple, Annotated, Union, TYPE_CHECKING # experimenting with typing
from numpy.typing import NDArray
from util.fastjet.jetfinderbase import JetFinderBase
from util.qol_utils.progress_bar import printProgressBarColor
from util.buffer.input import RootTreeLoader
from util.buffer.output import RootOutputBuffer
from util.misc.timing import profile_method, profile_block

import util.reconstruction.post_processing.utils.ghost_association as ghost_assoc
import util.reconstruction.post_processing.utils.softdrop as softdrop
import util.reconstruction.post_processing.utils.jhtagger as jhtagger
import util.reconstruction.post_processing.utils.jet_filter as jet_filter
import util.reconstruction.post_processing.utils.containment as containment
import util.reconstruction.post_processing.utils.simple_btag as simple_btag
import util.reconstruction.post_processing.utils.jet_energy_scale as jes

if(TYPE_CHECKING):
    from util.metadata.meta import MetaDataHandler

class JetFinder(JetFinderBase):
    """
    This class uses the Fastjet library to perform jet clustering, on some
    (arbitrary) set of inputs representing four-momenta of some objects.
    """

    def __init__(self, input_collections:List[str]=['StableTruthParticles'], jet_algorithm:str='anti_kt',radius:float=0.4, jet_name:str='AK04Jets', n_jets_max:int=-1,save_constituents:bool=True, fastjet_dir:Optional[str]=None,verbose:bool=False):

        super(JetFinder,self).__init__(fastjet_dir)
        self.status = False

        self.SetInputCollections(input_collections)
        self.jet_algorithm_name = jet_algorithm
        self.jet_name = jet_name
        self.radius = radius
        self.constituents_flag = save_constituents
        self.n_jets_max = n_jets_max # max number of jets to save per event (will be pt-ordered)
        self.n_constituents_max = 200 # max number of constituents to save per jet # TODO: Make configurable

        #TODO: Consider removing the "single_jet" functionality; it might make code maintenance harder?
        #      In principle one could do a pretty simple post-processing to whittle off the extra dimension if needed.
        self.single_jet = False # if true, self.n_jets_max = 1 & will remove the "number of jets" dimension (dim 1). Accessed by certain post-processors.
        if(self.n_jets_max == 1):
            self.single_jet = True # turn on if there's only 1 jet saved per event -- no real need for the extra dimension then

        self.fastjet_dir = fastjet_dir

        # Input buffer
        self.input_buffer = None

        # Output buffer
        self.output_buffer = RootOutputBuffer() # TODO: Make buffer size configurable. Larger sizes use more memory, but may be faster since we do fewer flushes and thus less I/O (depends on how good the flushing code is, shouldn't be open/closing files repeatedly!)

        self.input_collection_arrays = {}
        self.input_collection_arrays_cyl = {}
        self.input_collection_arrays_rapidity = {}
        self.constituent_indices_dict = None

        self.SetVerbosity(verbose)

        self.print_prefix = '\n\tJetFinder'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tRunning JetFinder:'
        self.progress_bar_suffix = 'Complete'

        self.error = False

        self._i = 0
        self.processors = [] # supposedly this is an example of an "observer pattern"
        self.processors_dict = {}

        self.metadata_handler = None

        # Stuff for dealing with I/O

        self.ntuple_file = None # The input HDF5 file -- also where output will ultimately be copied.
        self.output_file_tmp = None # The temporary output file -- write to this, then it'll be merged with input file. Allows having input and output files simultaneously read & modified.

        # Generate info on citations for algorithms
        self.citations = {}
        self._generate_citations()

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return

    def SetMetadataHandler(self,handler:'MetaDataHandler'):
        self.metadata_handler = handler

    def SetVerbosity(self,flag:bool):
        self.verbose = flag

    def SetNtupleFilename(self,file:str):
        """
        Sets the input n-tuple file.
        Also creates a name for the
        (temporary) output file.
        """
        self.ntuple_file = file

        self.output_file_tmp = self.ntuple_file
        extension = self.ntuple_file.split('.')[-1]
        self.output_file_tmp = '.'.join(self.ntuple_file.split('.')[:-1]) + '_postproc_tmp.' + extension

    def SetInputCollections(self,collections:List[str]):
        """
        This function prepares the actual key names that will be used
        to access collections' four-momenta from the (HDF5) n-tuple.

        We access both the Cartesian (Pmu) and cylindrical (Pmu_cyl),
        the former for constructing fastjet.PseudoJet objects, and the
        latter for accessing (eta,phi) which we can give to the PseudoJets
        to speed up certain operations where they would normally have to
        compute rapidity and phi (we can only use this *if* we're assuming
        eta ~= rapidity for the clustering inputs).
        """
        if(type(collections) != list):
            collections = [collections]

        self.input_collection_names = collections
        self.input_collection_names_Pmu = {collection: '{}.Pmu'.format(collection) for collection in collections}
        self.input_collection_names_Pmu_cyl = {collection: '{}.Pmu_cyl'.format(collection) for collection in collections}


        #NOTE: the keys for rapidity branches will be set later on

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

    def GetCitations(self):
        return self.citations

    # def _input_consistency_check(self):
    #     f = h5.File(self.ntuple_file,'r')
    #     keys = list(f.keys())
    #     f.close()

    #     cleaned_collections = []
    #     for collection in self.input_collection_names_Pmu:
    #         # we will be using the Cartesian versions of each collection for clustering
    #         found = collection in keys
    #         if(not found):
    #             self._print('Warning: Did not find key {} in file {}. Disabling as input...'.format(collection,self.ntuple_file))
    #         else:
    #             cleaned_collections.append(collection)
    #     self.input_collection_names_Pmu = cleaned_collections
    #     if(len(self.input_collection_names_Pmu)==0):
    #         self._print('Error: No input collections.')
    #         return False
    #     return True

    def _generate_citations(self):
        """
        Fills in citations (in BibTex format) for FastJet.
        Additional post-processors can add to this.
        """
        key = 'FastJet'
        if(key not in self.citations.keys()):
            self.citations[key] = [
                """
@article{Cacciari:2011ma,
    author = "Cacciari, Matteo and Salam, Gavin P. and Soyez, Gregory",
    title = "{FastJet User Manual}",
    eprint = "1111.6097",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "CERN-PH-TH-2011-297",
    doi = "10.1140/epjc/s10052-012-1896-2",
    journal = "Eur. Phys. J. C",
    volume = "72",
    pages = "1896",
    year = "2012"
}
                """,
                """
@article{Cacciari:2005hq,
    author = "Cacciari, Matteo and Salam, Gavin P.",
    title = "{Dispelling the $N^{3}$ myth for the $k_t$ jet-finder}",
    eprint = "hep-ph/0512210",
    archivePrefix = "arXiv",
    reportNumber = "LPTHE-05-32",
    doi = "10.1016/j.physletb.2006.08.037",
    journal = "Phys. Lett. B",
    volume = "641",
    pages = "57--61",
    year = "2006"
}
                """
            ]

        # Now, fetch citations for any existing post-processors.
        for post_proc in self.processors:
            for k,v in post_proc.GetCitations().items():
                if(k) not in self.citations.keys():
                    self.citations[k] = v

    def Initialize(self):
        """
        Reads in the input HDF5 file, and places the required arrays in memory.
        """
        # TODO: May want to consider chunking things and using a buffer?
        # Memory usage will scale better for larger files.
        # However, need to be careful if input and output are the same file,
        # since writing maybe needs to keep it open in order to avoid lots of
        # opening/closing that will slow down the program.

        # If already initialized, no need to do it again.
        if(self.status):
            return

        # # Input consistency check.
        # if(not self._input_consistency_check()):
        #     self.status = False
        #     return

        # Initialize fastjet.
        self._initialize_fastjet()
        if(not self.fastjet_init_flag):
            self.status = False
            return

        # Initialize the fastjet jet definition
        self._initialize_jet_definition()

        # Read in the input 4-momenta from the input file.

        self._fetch_inputs()

        # Also fetch rapidity & phi, for potentially speeding up some FastJet computations.
        self._fetch_rapidity() # TODO: Fix this

        # Get maximum size of jet inputs.
        n_max = self._get_max_input_size()
        self._initialize_pseudojets(n_max)

        # Optional modification of initialize. May be harnessed by some special configurations.
        self._modifyInitialization()

        # With nevents defined, we can initialize the buffer.
        self._initializeBuffer()

        # Now (re)generate citations, will pull in any additions from post-processors that have been added on.
        self._generate_citations()

        self.status = True

    def _fetch_inputs(self):

        self.input_buffer = RootTreeLoader(self.ntuple_file,'hepdata4ml_tree') #TODO: Dynamic tree name?
        self.input_buffer.load()
        for collection_name,key in self.input_collection_names_Pmu.items():
            self.input_buffer.read_branch(key)
        for collection_name,key in self.input_collection_names_Pmu_cyl.items():
            self.input_buffer.read_branch(key)

        self.nevents = self.input_buffer.t.GetEntries()

        return

    def _fetch_rapidity(self):
        """
        Fetches the rapidity of the jet clustering inputs.
        This can be explicitly passed on to FastJet, and should speed
        up the clustering -- which we ought to do if we've already
        spent time computing it.

        Note: We may use eta instead of rapidity and thus implicitly assume
        the input 4-vecs to be massless (as they often are). However,
        we try to fetch any existing "Rapidity"/"Rap"/"Y" branch first,
        in case it exists.
        """

        # For rapidity, fetch rapidity or pseudorapidity based on what is available.
        # Thus we will store keys in "self.input_collection_names_rapidity", where the actual rapidity keys
        # are mapped to contents of self.input_collection_names.
        # Note that self.input_collection_names_Pmu(_cyl) etc are built from self.input_collection_names.
        full_keys = list(self.input_buffer.keys) # gives all available branch names in the input TTree
        self.input_collection_names_rapidity = {}
        for key in self.input_collection_names:
            potential_keys = ['{}.{}'.format(key,x) for x in ['Rapidity','Rap','Y','Pmu_cyl']] # last is the fall-back
            for key2 in potential_keys:
                if(key2 in full_keys):
                    self.input_collection_names_rapidity[key] = key2
                    break
        assert(len(self.input_collection_names_rapidity.keys()) == len(self.input_collection_names))

        # Now, ensure that the necessary collections are loaded
        for key,key2 in self.input_collection_names_rapidity.items():
            self.input_buffer.read_branch(key2)
        return

    def _get_max_input_size(self):
        """
        Determines the maximum number of pseudojets we'll need.
        Note that this likely overshoots the maximum, because
        it sums up the maximum of each collection; it's not
        guaranteed that this maxima all come from the same event.
        """
        # Access the TTree directly from self.input_buffer
        sizes = {}
        for key in self.input_collection_names:
            branch_name = '{}.N'.format(key)
            sizes[key] = self.input_buffer.t.GetMaximum(branch_name)

        return int(np.sum(np.stack(list(sizes.values()))))

    def _load_data(self,event_index=None):
        """
        This function fills self.input_collection_arrays etc.
        """
        if(event_index is None):
            event_index = self._i

        self.input_buffer.set_entry(event_index)

        #TODO: Likely have to fix some things with the ghost associator
        for collection_name,key in self.input_collection_names_Pmu.items(): # NOTE: Using self.input_collections_array.keys() can be dangerous, due to modifications/additions to keys by things like GhostAssociation(). Those should not touch self.input_collections, for this reason.
            self.input_collection_arrays[collection_name] = self.input_buffer[key]

        for collection_name,key in self.input_collection_names_Pmu_cyl.items(): # NOTE: Using self.input_collections_array.keys() can be dangerous, due to modifications/additions to keys by things like GhostAssociation(). Those should not touch self.input_collections, for this reason.
            self.input_collection_arrays_cyl[collection_name] = self.input_buffer[key]

        # Also take care of rapidity. Here, we have to keep in mind the case where we've pointed
        # the rapidity collection at a cylindrical Pmu key (Pmu_cyl), in which case we're
        # linked to a branch carrying not the rapidity but the whole four-momentum; in that
        # case we need to peel off the eta component (it'll be pseudorapidity in this case).
        #
        # We'll actually *handle* this complication in _set_inputs(), however.
        for collection_name,key in self.input_collection_names_rapidity.items():
            self.input_collection_arrays_rapidity[collection_name] = self.input_buffer[key]

        # print('Loaded event {}'.format(event_index))
        # key = self.input_collection_names_Pmu_cyl[0]
        # print('\tPrinting buffer for {}'.format(key))
        # print(self.input_collection_arrays_cyl[key])

        return

    def _set_inputs(self):

        # We want to vstack the input_collection_arrays, but have to consider the edge
        # case where one of them is empty, in which case it'll be "{}". This will cause
        # dimensionality issues with vstack if we do things naively.
        self.SetInputs(np.vstack([self.input_collection_arrays[cname] for cname in self.input_collection_names_Pmu.keys() if len(self.input_collection_arrays[cname]) > 0])) # NOTE: Using self.input_collections_array.keys() can be dangerous, due to modifications/additions to keys by things like GhostAssociation(). Those should not touch self.input_collections, for this reason.
        self.SetInputsCylindrical(np.vstack([self.input_collection_arrays_cyl[cname] for cname in self.input_collection_names_Pmu_cyl.keys() if len(self.input_collection_arrays_cyl[cname]) > 0]))

        # For rapidity, there is the complication that we might be reading a "Pmu_cyl" branch, in case
        # we really just want its pseudorapidity component: slicing like [:,1] but these are cppyy.gbl.std.vector,
        # not numpy arrays, so we have to do it correctly.

        # Collect views/slices
        with profile_block('rapidity setting'):
            rapidity_views = []
            for cname in self.input_collection_names_rapidity.keys():
                coll = self.input_collection_arrays_rapidity[cname]
                rapidity = np.asarray(coll)
                rapidity_views.append(rapidity[:, 1] if rapidity.ndim > 1 else rapidity)
            self.SetRapidity(np.concatenate(rapidity_views, axis=0))

    @profile_method('JetFinder.Process')
    def Process(self):
        self.Initialize()

        if(self.verbose):
            self._print('Process: Number of events = {}'.format(self.nevents))
            printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        # Event loop
        for self._i in range(self.nevents):

            # Clear the user info, which marks jets in this event.
            # May be harnessed by some special configurations.
            self.ClearUserInfo()

            self._load_data() # load this event, fill

            # Gather the different input collections together, into one array of four-momenta.
            self._set_inputs()

            # Optional modification of inputs. May be harnessed by some special configurations.
            self._modifyInputs()

            self._clusterJets() # fills self.jets_dict

            # optionally extract information on jet constituents
            if(self.constituents_flag):
                self._fetchJetConstituents() # fills self.constituent_vectors, self.constituent_indices

            # Now run post-processing methods -- incl. writing to buffer.
            # These may modify the jets themselves, so we only call _writeToBuffer() after.
            # (If these methods need to access the jets, they will access them directly and not from buffer).
            self.PostProcess()

            # Pt-sort the jets, and truncate to fixed length given by self.n_jets_max
            self._ptSort(truncate=True)

            # # Pass the pt ordering to the buffer; will be used for writing out any branches
            # # that post-processors have currently stored in self.output_buffer.buffer_dicts.
            # # (This is the preferred way for post-processors to write out jet- and constituent-level
            # #  branches, because the final ordering of the jets isn't known when they're writing
            # #  to buffer).
            # self._print('_i = {}, setting ordering = {}'.format(self._i,self.jet_ordering))
            # self.output_buffer.set_ordering(self.jet_ordering)

            # now write jets to buffer
            self._writeToBuffer()

            if(self.verbose):
                printProgressBarColor(self._i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)
        return

    def _modifyInitialization(self):
        """
        Initialize the post-processors.
        """
        for processor in self.processors:
            processor.ModifyInitialization(self)
        return

    def _modifyInputs(self):
        for processor in self.processors:
            processor.ModifyInputs(self)
        return

    def PostProcess(self):
        """
        Sequentially run the post-processors, to modify jets (and possibly constituents).
        This may modify self.jets_dict and its derived quantities, and/or create new
        branches in the output buffer.
        """
        # NOTE: Used to "interleave" post-processors with main jet clustering method,
        #       but doing this fully sequentially is probably better; allows for nicer
        #       and clearer interplay between the post-processors (each can access outputs
        #       of previous ones fully).
        #
        #       Note that we *do* still effectively interleave the ModifyInitialization()
        #       and ModifyInputs() methods with the jet clustering.
        for processor in self.processors:
            processor.ModifyJets(self)
            processor.ModifyConstituents(self)
            processor.ModifyWrite(self)
        return

    def Flush(self, output_file=None):
        """
        This function simply finishes the writing of our data buffer
        to the temporary output file, by doing a final flush.
        Then, the temporary output file is merged into the final output
        file (which typically is the input file -- we're just appending to it).
        It also writes some metadata, which will ultimately propagate
        to the output file.
        """
        # Final flush of the buffer.
        self.output_buffer.flush()

        # Now, handle the final output file.
        if(output_file is None):
            output_file = self.ntuple_file
        self.output_buffer.close(output_file) # <- merges the buffer's temporary output into ouput_file, deletes tmp output

        self._writeMetadata()
        return

    def _writeMetadata(self):
        """
        Writes metadata to the output file.
        """

        if(self.metadata_handler is None):
            self._print('Unable to write metadata; no handler was provided.')
            return

        metadata = self.metadata_handler.GetMetaData()

        # Add information on jet input collections.
        key = 'Metadata.JetCollections.InputCollections'
        if(key not in metadata.keys()):
            metadata[key] = {}
        metadata[key][self.jet_name] = [x for x in self.input_collection_names]

        # Also add metadata on citations for algorithms. This is stored as a list of strings; it is not separated by jet_name,
        # as that level of granularity is probably not useful.
        key = 'Metadata.Citations'
        self.metadata_handler.AddCitations(self.GetCitations())

    def __call__(self,hepmc_file:str,h5_file:str,output_file:Optional[str]=None,verbose:Optional[bool]=None, copts:int=9, key:Optional[str]=None):
        """
        Note that this uses a generic signature, although some arguments may be unused.
        This is to keep the structure similar between different post-processors.
        """
        # should consider making the various post-processors inherit from a single parent class!

        if(verbose is not None): self.SetVerbosity(verbose)
        self.SetNtupleFilename(h5_file)

        self.Process()
        self.Flush(output_file)

        return self.ntuple_file

    def _initializeBuffer(self):

        self.output_buffer.SetFilename(self.output_file_tmp)
        self.output_buffer.SetCloneTree(self.input_buffer.GetTree()) # clone the input n-tuple TTree structure -- will be filled as we loop through and call TTree::Fill()

        dim0 = 1
        dim1 = 2
        dim2 = 3

        if(self.single_jet): # eliminate the "number of jets" dimension
            dim0 = 0
            dim1 = 1
            dim2 = 2

        self.output_buffer.create_array('{}.N'.format(self.jet_name),dtype=np.dtype('i4'))
        self.output_buffer.create_array('{}.Pmu'.format(self.jet_name),ndim=dim1,dtype=np.dtype('f8'))
        self.output_buffer.create_array('{}.Pmu_cyl'.format(self.jet_name),ndim=dim1,dtype=np.dtype('f8'))
        if(self.constituents_flag):
            self.output_buffer.create_array('{}.Constituents.N'.format(self.jet_name),ndim=dim0,dtype=np.dtype('i4'))
            self.output_buffer.create_array('{}.Constituents.Pmu'.format(self.jet_name),ndim=dim2,dtype=np.dtype('f8'))
            self.output_buffer.create_array('{}.Constituents.Pmu_cyl'.format(self.jet_name),ndim=dim2,dtype=np.dtype('f8'))

            # Also create buffers corresponding to jet constituents' indices w.r.t. the collections they were pulled from.
            # Note that a jet may have used multiple collections -- so we'll keep track of the index of the collection that
            # a constituent came from, as well as its index *within* that collection.
            self.output_buffer.create_array('{}.Constituents.Collection'.format(self.jet_name),ndim=dim1,dtype=np.dtype('i4'))
            self.output_buffer.create_array('{}.Constituents.Collection.Index'.format(self.jet_name),ndim=dim1,dtype=np.dtype('i4'))

        return

    @profile_method('JetFinder._computeConstituentIndices')
    def _computeConstituentIndices(self):
        # Precompute collection boundaries once
        n_per_collection = [len(self.input_collection_arrays[cname]) for cname in self.input_collection_names_Pmu.keys()]
        cumulative_lengths = np.cumsum([0] + n_per_collection)

        self.constituent_indices_dict = {}
        for i, jet in self.jets_dict.items():
            raw_indices = np.array([pj.user_index() for pj in jet.constituents()])
            # Vectorized conversion for all indices at once
            collection_indices = np.searchsorted(cumulative_lengths[1:], raw_indices, side='right')
            local_indices = raw_indices - cumulative_lengths[collection_indices]

            # Combine into pairs
            constituent_indices = np.array(list(zip(collection_indices, local_indices)),dtype=np.dtype('i4'))

            # Store the results
            self.constituent_indices_dict[i] = constituent_indices
        return

    @profile_method('JetFinder._writeToBuffer')
    def _writeToBuffer(self,event_index:Optional[int]=None):

        if(event_index is None):
            event_index = self._i

        if(len(self.jets_dict) == 0):
            self.output_buffer.flush() # should write an empty entry
            return
        njets = len(self.jet_vectors)

        # Fill jet information in the buffer.
        self.output_buffer.set('{}.N'.format(self.jet_name),event_index,njets)

        self._load_data()

        # TODO: Maybe later clean this up a bit? Have to deal with special case of "single_jet = True".
        if(self.single_jet):
            idx = self.jet_ordering[0]
            self.output_buffer.set('{}.Pmu'.format(self.jet_name),event_index,self.jet_vectors[idx])
            self.output_buffer.set('{}.Pmu_cyl'.format(self.jet_name),event_index,self.jet_vectors_cyl[idx])

        else:
            self.output_buffer.set('{}.Pmu'.format(self.jet_name),event_index,self.jet_vectors)
            self.output_buffer.set('{}.Pmu_cyl'.format(self.jet_name),event_index,self.jet_vectors_cyl)

        # Fill the jet constituent information.
        if(self.constituents_flag):

             # Figure out the collections and indices of the constituents
            self._computeConstituentIndices()

            if(self.single_jet):
                idx = self.jet_ordering[0]
                self.output_buffer.set('{}.Constituents.N'.format(self.jet_name),event_index,len(self.constituent_vectors[idx]))

                self.output_buffer.set('{}.Constituents.Pmu'.format(self.jet_name),event_index,self.constituent_vectors[idx])
                self.output_buffer.set('{}.Constituents.Pmu_cyl'.format(self.jet_name),event_index,self.constituent_vectors_cyl[idx])
                self.output_buffer.set('{}.Constituents.Collection'.format(self.jet_name),event_index,self.constituent_indices_dict[idx][:,0])
                self.output_buffer.set('{}.Constituents.Collection.Index'.format(self.jet_name),event_index,self.constituent_indices_dict[idx][:,1])

            else:
                self.output_buffer.set('{}.Constituents.N'.format(self.jet_name),event_index,{i:len(self.constituent_vectors[i]) for i in self.constituent_vectors.keys()})
                self.output_buffer.set('{}.Constituents.Pmu'.format(self.jet_name),event_index,self.constituent_vectors)
                self.output_buffer.set('{}.Constituents.Pmu_cyl'.format(self.jet_name),event_index,self.constituent_vectors_cyl)
                self.output_buffer.set('{}.Constituents.Collection'.format(self.jet_name),event_index,{key:val[:,0] for key,val in self.constituent_indices_dict.items()})
                self.output_buffer.set('{}.Constituents.Collection.Index'.format(self.jet_name),event_index,{key:val[:,1] for key,val in self.constituent_indices_dict.items()})

        # Lastly, we tell the buffer what is the jet ordering for this event.
        # Next time it flushes (at the top of this loop on the next iteration, or on explicit flush),
        # it will use this to determine how to actually sort the dictionaries of jets/constituents we've provided.
        self.output_buffer.set_ordering(self.jet_ordering)

        return

    # NOTE: Will define various functions for performing some modifications to clustering or post-processing of results.
    #       Would be nice to put this somewhere else, but I also want to keep the configuration simple.
    #       These will be member functions that return self, so you can do "constructor().function()" instead of just "constructor()"
    #       and in this way chain together a complex configuration without the constructor having to take a huge number of args.

    def PtFilter(self,pt_min:Annotated[float,"GeV"]=15.):
        self.processors.append(jet_filter.PtFilter(pt_min))
        return self

    def EtaFilter(self,eta_max:float=2.):
        self.processors.append(jet_filter.EtaFilter(eta_max))
        return self

    def Leading(self):
        self.processors.append(jet_filter.Leading())
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

    def Softdrop(self,z_cut:float,beta:float):
        """
        This function performs the softdrop algorithm.

        Returns self, so this can be chained with the constructor.
        """
        self.processors.append(softdrop.Softdrop(z_cut,beta))
        return self

    def IteratedSoftdrop(self,z_cut:float,beta:float,dR_cut:float,max_depth:int=10,mode:str='tag'):
        """
        This function performs the iterated softdrop algorithm.

        Returns self, so this can be chained with the constructor.
        """
        self.processors.append(softdrop.IteratedSoftdrop(z_cut,beta,dR_cut,max_depth,mode))
        return self

    def JohnsHopkinsTagger(self,delta_p:float=0.1,delta_r:float=0.19,cos_theta_W_max:float=0.7,top_mass_range:Annotated[Tuple[float,float],"GeV"]=(150.,200.),W_mass_range:Annotated[Tuple[float,float],"GeV"]=(65.,95.), mode:str='filter',tag_name:Optional[str]=None):
        """
        This function performs top-tagging via
        the Johns Hopkins top tagger.

        With mode=='filter', it will filter out non-associated jets.
        With mode=='tag', it will save a ghost association flag to the output.

        Returns self, so this can be chained with the constructor.
        """

        self.processors.append(jhtagger.JohnsHopkinsTagger(delta_p,delta_r,cos_theta_W_max,top_mass_range,W_mass_range,mode,tag_name))
        return self

    def Containment(self,truth_key:str,truth_indices:Union[int,List[int],NDArray],delta_r:Optional[float]=None,use_rapidity:bool=True, mode:str='tag',tag_name:Optional[str]=None):
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

    def TrackCountingBTag(self,track_key:str,track_pt_min:Annotated[float,"GeV"]=1., delta_r:float=0.3, track_ip_max:Annotated[float,"mm"]=2., sig_min:float=6.5, ntracks:int=3, use_3d:bool=False, mode:str='tag',tag_name:Optional[str]=None):
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

    def Leading(self):
        self.processors.append(jet_filter.Leading())
        return self

    def JESCalibration(self,formula:str):
        """
        This function performs a jet energy scale calibration of the jets.
        It directly modifies the jets, so their momenta (and constituent momenta)
        are calibrated. Also writes the per-jet calibration factor to a new branch,
        with which one can determine the pre-calibrated momenta.

        Returns self, so this can be chained with the constructor.
        """
        self.processors.append(jes.JetEnergyScale(formula))
        return self

    def JESCalibrationATLAS(self):
        """
        This function performs a jet energy scale calibration of the jets,
        corresponding to the one in the ATLAS Delphes card.
        Returns self, so this can be chained with the constructor.
        """
        self.processors.append(jes.JetEnergyScale('atlas'))
        return self

    def JESCalibrationCMS(self):
        """
        This function performs a jet energy scale calibration of the jets,
        corresponding to the one in the default CMS Delphes card.
        Returns self, so this can be chained with the constructor.
        """
        self.processors.append(jes.JetEnergyScale('cms'))
        return self

class TruthJetFinder(JetFinderBase):
    """
    A simple jet-finding class, for use with event filters.
    """

    def __init__(self, jet_algorithm:str='anti_kt',radius:float=0.4, jet_name:str='AK04Jets', n_jets_max:Optional[int]=None,fastjet_dir:Optional[str]=None):

        super(TruthJetFinder,self).__init__(fastjet_dir)

        self.status = False

        self.jet_algorithm_name = jet_algorithm
        self.jet_name = jet_name
        self.radius = radius
        self.n_jets_max = n_jets_max # max number of jets to save per event (will be pt-ordered)
        self.n_constituents_max = 200 # max number of constituents to save per jet

        self.fastjet_dir = fastjet_dir
        self.fastjet_init_flag = False

        self.print_prefix = '\n\tTruthJetFinder'

        self.error = False

        self._i = 0

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def Initialize(self):

        # If already initialized, no need to do it again.
        if(self.status):
            return

        # Initialize fastjet.
        self._initialize_fastjet()
        if(not self.fastjet_init_flag):
            self.status = False
            return

        # Initialize the fastjet jet definition
        self._initialize_jet_definition()

        self.status = True

    def Process(self,input_vecs):
        self.input_vecs = input_vecs
        self._clusterJets()
        self._ptSort()

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return
