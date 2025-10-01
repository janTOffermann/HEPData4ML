import sys, operator
import numpy as np
import ROOT as rt
from util.fastjet.setup import FastJetSetup, FastJetInterfaceSetup
from typing import Optional, Any
from util.misc.timing import profile_method, profile_block

class JetFinderBase:
    """
    This is a base class, for using the Fastjet library to perform jet clustering.
    It is leveraged by the JetFinder class in utils/post_processing/jets.py .
    This base class exists to avoid some possible import loops.
    It lacks some useful functionality that JetFinder includes.
    """

    def __init__(self,fastjet_dir:Optional[str]=None):

        self.jet_algorithm_name = ''
        self.jet_name = ''
        self.radius = 0.4

        self.print_prefix = '\n\tJetFinderBase'

        self.fastjet_dir = fastjet_dir
        self.fastjet_init_flag = False

        # some fastjet-specific vars, for internal usage
        self.jet_algorithm = None
        self.jetdef = None
        self.cluster_sequence = None
        # self.jets = None
        self.jets_dict = None
        self.jet_vectors = None
        self.jet_vectors_cyl = None
        self.constituent_vectors = None
        self.constituent_vectors_cyl = None
        self.constituent_indices = None
        self.user_info = None
        self.pt_sorting = None # allows access to the sorting array
        self.n_jets_max = None
        self.n_constituents_max = None

        self.input_vecs = None
        self.input_vecs_cyl = None
        self.rapidity = None
        self.configurator = None

        self.jet_ordering = None

        # Buffer containing FastJet::PseudoJet objects -- better to use
        self.pseudojets = None
        self.pseudojet_init_flag = False

    def _initialize_pseudojets(self, size=1000, force=False):
        if(self.pseudojet_init_flag and not force):
            return

        self._initialize_fastjet()
        import fastjet as fj
        self.pseudojets = [fj.PseudoJet() for i in range(size)]

        self.pseudojet_init_flag = True
        return

    def _initialize_fastjet(self):

        if(self.fastjet_init_flag):
            return

        if(self.fastjet_dir is None):
            if(self.configurator is None):
                self._print('Error: Fastjet being requested, but this JetFinderBase has no self.fastjet_dir nor self.configurator . Fastjet import will not work.')
                return # bad
            else: # Fetch fastjet directory from configurator. This is foreseen as the "typical" usage.
                self.fastjet_dir = self.configurator.GetFastjetDirectory()

        # Now, if possible we initialize fastjet
        # if(self.fastjet_dir is not None):
        self._setupFastJet()
        sys.path.append(self.fastjet_dir)
        import fastjet as fj # This is where fastjet is really imported & cached by Python, but there may be other import statements peppered throughout since this has limited scope.
        self.fastjet_init_flag = True
        return

    def SetNConstituentsMax(self,n):
        self.n_constituents_max = n

    def SetRadius(self,radius):
        self.radius = radius

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def SetInputs(self,vecs):
        self.input_vecs = vecs

    def SetInputsCylindrical(self,vecs):
        self.input_vecs_cyl = vecs

    def SetRapidity(self,val):
        self.rapidity = val

    def _setupFastJet(self):
        verbose = self.configurator.GetPrintFastjet()
        self.fastjet_setup = FastJetSetup(self.configurator.GetFastjetDirectory(),full_setup=True,verbose=verbose)
        self.configurator.SetPrintFastjet(False)
        self.fastjet_dir = self.fastjet_setup.GetPythonDirectory()
        return

    def _parse_jet_algorithm(self):
        self._initialize_fastjet()
        import fastjet as fj

        self.jet_algorithm = None
        for key in ['anti_kt','anti kt','anti-kt']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = fj.antikt_algorithm
                return True

        for key in ['kt']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = fj.kt_algorithm
                return True

        for key in ['c/a','cambridge','aachen']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = fj.cambridge_aachen_algorithm
                return True

        return False

    def _initialize_jet_definition(self):
        self._initialize_fastjet()
        import fastjet as fj
        # Print the FastJet banner -- it's unavoidable (other packages don't do this!).
        fj.ClusterSequence.print_banner()

        # determine what jet algorithm to use
        self._parse_jet_algorithm()

        self.jetdef = fj.JetDefinition(self.jet_algorithm, self.radius)
        return

    @profile_method('JetFinderBase._clusterJets')
    def _clusterJets(self):
        self._initialize_fastjet()
        import fastjet as fj # hacky, but will work because _setupFastJet() was run in __init__()

        # Quick check to make sure our PseudoJet buffer is large enough.
        # Ideally, we should prepare this upstream to avoid having to re-initialize.
        n_pseudojets = len(self.input_vecs)
        if(n_pseudojets > len(self.pseudojets)):
            self._initialize_pseudojets(n_pseudojets,force=True) # TODO: Consider adding a safety factor?

        # Now set the PseudoJet momenta and indices -- the latter for tracing them through jet clustering.
        # vecs has format (E,px,py,pz) -- FastJet uses (px,py,pz,E) so we must rearrange. Faster than np.roll.
        has_rapidity = (self.rapidity is not None)
        for i,x in enumerate(self.input_vecs):
            self.pseudojets[i].reset(x[1],x[2],x[3],x[0]) # NOTE: This should reset indices -- and user info (?).
            self.pseudojets[i].set_user_index(i)

            # If we've supplied (rapidity,phi) of inputs to JetFinderBase, we can pass these to FastJet to avoid
            # having it recompute these quantities internally.
            if(has_rapidity):
                self.pseudojets[i].set_cached_rap_phi(self.rapidity[i],self.input_vecs_cyl[i,2])

        # Attach any optional information to the pseudojet objects. This can be leveraged by other classes
        # or extensions.
        if(self.user_info is not None):
            for idx, val in self.user_info.items(): # in practice, val will be a dictionary itself -- allows for attaching multiple things
                self.pseudojets[idx].set_python_info(val)

        with profile_block('JetFinderBase._clusterJets - ClusterSequence'): # Useful for profiling -- this time block is largely non-negotiable.
            self.cluster_sequence = fj.ClusterSequence(self.pseudojets[:n_pseudojets], self.jetdef) # member of class, otherwise goes out-of-scope when ref'd later
        self.jets_dict = {i:jet for i,jet in enumerate(self.cluster_sequence.inclusive_jets())} # NOTE: Repeated calls to ClusterSequence::inclusive_jets() seems OK, I think it is just an accessor.
        self.jet_ordering = np.arange(len(self.jets_dict))
        self.pt_sorting = self.jet_ordering
        self._jetsToVectors()
        return

    def _jetsToVectors(self):
        """
        This function fills self.jet_vectors and self.jet_vectors_cyl, to contain the four-momenta
        of whatever jets are currently in self.jets.
        """
        self.jet_vectors = {i:np.array([jet.e(), jet.px(), jet.py(),jet.pz()]) for i,jet in self.jets_dict.items()}
        self.jet_vectors_cyl = {i:np.array([jet.pt(), jet.eta(), jet.phi(),jet.m()]) for i,jet in self.jets_dict.items()}

    @profile_method('JetFinderBase._ptSort')
    def _ptSort(self, truncate=False):
        """
        Sorts jets by decreasing pT, and truncates to
        take only the first self.n_jets_max jets.

        This is accomplished by modifying self.jet_ordering.
        """

        if(len(self.jets_dict) < 1): # no jets -> nothing to do (for case of 1 jet, we'll modify some stuff below)
            return
        elif(len(self.jets_dict) == 1): # 1 jet -> not much to do, just make sure self.jet_ordering reflects this
            self.jet_ordering = [self.jet_ordering[0]]
            return

        # For the jet pt, access self.jet_vectors_cyl instead of calling PseudoJet::pt()
        # on the contents of jets_dict; this should be faster as we avoid having FastJet
        # recompute these quantities.
        jet_pt = np.array([self.jet_vectors_cyl[i][0] for i in self.jet_ordering])
        is_sorted = np.all(jet_pt[:-1] >= jet_pt[1:])
        needs_truncation = (self.n_jets_max is not None and truncate and len(jet_pt) > self.n_jets_max)

        # Early return if already sorted and no truncation needed
        if is_sorted and not needs_truncation:
            # In principle, we could check that self.pt_sorting is not None,
            # but calling _clusterJets() will ensure this is filled.
            # Still need to set pt_sorting if it doesn't exist
            return

        with profile_block('JetFinderBase._ptSort - Sort'):
            if (not is_sorted):
                self.pt_sorting = np.argsort(-jet_pt)

            if(needs_truncation):
                self.pt_sorting = self.pt_sorting[:self.n_jets_max]

            if(len(self.pt_sorting) == 1):
                self.jet_ordering = [self.jet_ordering[self.pt_sorting[0]]]
            else:
                self.jet_ordering = list(operator.itemgetter(*self.pt_sorting)(self.jet_ordering))

            # If we've truncated, we'll need to update some dictionaries under-the-hood to reflect this.
            if(needs_truncation):
                # remove entries from jets_dict, that correspond with entries in jet_ordering that have been dropped
                self._updateJetDictionary()

                # We will also recompute the jet vectors, to account for any that have been dropped.
                # Note that due to the dictionary-based approach, we don't have to recompute this if
                # the jet ordering has simply changed.
                self._jetsToVectors()

                # Also refresh constituents. Again, this only needs to be called if jets were dropped,
                # since its a dictionary so a simple reordering of the jets in self.jet_ordering does
                # not necessitate any change.
                self._fetchJetConstituents()
            return

    def _updateJetDictionary(self):
        """
        To be used when self.jet_ordering is updated.
        """
        for key in list(self.jets_dict.keys()):
            if(key not in self.jet_ordering):
                del self.jets_dict[key]
        return

    @profile_method('JetFinderBase._fetchJetConstituents')
    def _fetchJetConstituents(self):
        results = {i:self._fetchJetConstituentsSingle(i, jet, self.n_constituents_max) for i,jet in self.jets_dict.items()}
        self.constituent_vectors = {i:self.input_vecs[x] for i,x in results.items()}
        self.constituent_vectors_cyl = {i:self.input_vecs_cyl[x] for i,x in results.items()}
        self.constituent_indices = {i:x for i,x in results.items()}

    def _fetchJetConstituentsSingle(self, key, jet, n_constituents=-1):
        """
        Returns indices of the jet constituents, w.r.t. the
        self.input_vecs list. The indices have been pt-sorted
        and truncated if requested.
        """
        if not jet.has_constituents():
            return np.empty((0, 4)), np.empty((0, 4)), np.empty((0, 4))

        constituents = jet.constituents()
        n = len(constituents)

        if n_constituents is not None and n_constituents > 0:
            max_constituents = min(n_constituents, n)
        else:
            max_constituents = n

        # Get the user indices of the constituents.
        # The input 4-vectors were indexed sequentially, so we can
        # use this to look them up in our original inputs, thus avoiding
        # calls to Fastjet::Pseudojet. Especially useful for the coordinates
        # that it internally recalculates -- pt, eta, phi and m (in case we've
        # supplied rap/eta/phi, those might not be recalculated anyway).
        indices = np.array([constituent.user_index() for constituent in constituents],dtype=int)

        # Sort on pt, truncate if necessary.
        pt = self.input_vecs_cyl[indices,0]

        # For the indices, we can directly return the sorted indices
        # since, by construction, the "indices" are just a 0-indexed range.
        result = np.argsort(-pt)[:max_constituents]

        # print('######')
        # print('# {}'.format(key))
        # print('Constituent indices: ', result)
        # print('pt: ',pt[result])
        # print('######')
        # print()

        return result

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return

class ParticleInfo(object):
    """Illustrative class for use in assigning pythonic user information
    to a PseudoJet.
    """
    def __init__(self, particle_index, status, pdg_id=0):
        self.particle_index = particle_index
        self.status = status
        self.pdg_id = pdg_id

    def __str__(self):
        return "particle_index={0}, status={1}, pdg_id={2}".format(
            self.particle_index, self.status, self.pdg_id)






class FastJetFinderBase:
    """
    This is a base class, for using the Fastjet library to perform jet clustering.
    It is leveraged by the JetFinder class in utils/post_processing/jets.py .
    Based on the JetFinderBase class, but using a custom fastjet interface
    that might speed things up.
    """

    def __init__(self,fastjet_dir:Optional[str]=None):

        self.setup = FastJetInterfaceSetup(fastjet_dir)
        self.init = False
        self.fastjet_interface = None

        self.jet_algorithm_name = ''
        self.jet_name = ''
        self.radius = 0.4

        self.print_prefix = '\n\tFastJetFinderBase'

        self.fastjet_dir = fastjet_dir
        self.fastjet_init_flag = False

        # some fastjet-specific vars, for internal usage
        self.jet_algorithm = None
        self.jetdef = None
        self.cluster_sequence = None
        self.jet_vectors = None
        self.jet_vectors_cyl = None
        self.constituent_vectors = None
        self.constituent_vectors_cyl = None
        self.constituent_indices = None
        self.user_info = None
        self.pt_sorting = None # allows access to the sorting array
        self.n_jets_max = None
        self.n_constituents_max = None

        self.input_vecs = None
        self.input_vecs_cyl = None
        self.configurator = None

        self.jet_ordering = None

    def _initialize_fastjet(self):
        """
        Here, we make sure that our FastjetInterface
        library is built and can be loaded with PyROOT.
        """
        if(self.init):
            return
        self.setup.FullPreparation()
        self.init = True

        self.fastjet_interface = rt.FastjetInterface.JetFinder()
        return

    def SetNConstituentsMax(self,n):
        self.n_constituents_max = n

    def SetRadius(self,radius):
        self.radius = radius

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def SetInputs(self,vecs):
        self.input_vecs = np.asarray(vecs)

    def SetInputsCylindrical(self,vecs):
        self.input_vecs_cyl = np.asarray(vecs)

    def _parse_jet_algorithm(self):

        self.jet_algorithm = None
        for key in ['anti_kt','anti kt','anti-kt']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = 'antikt'
                return True

        for key in ['kt']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = 'kt'
                return True

        for key in ['c/a','cambridge','aachen']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = 'ca'
                return True

        return False

    def _initialize_jet_definition(self):
        self._initialize_fastjet()

        # # Print the FastJet banner -- it's unavoidable (other packages don't do this!).
        # fj.ClusterSequence.print_banner()

        # determine what jet algorithm to use
        self._parse_jet_algorithm()
        self.fastjet_interface.SetJetAlgorithm(self.jet_algorithm)
        self.fastjet_interface.SetRadius(self.radius)
        return

    @profile_method('FastJetFinderBase._clusterJets')
    def _clusterJets(self):
        self._initialize_fastjet()

        # Set the inputs. We have a few different ways of doing this.
        # The SetInputsWithEtaPhiM will take (eta,phi,m) in addition as a ways
        # of internally determining input vectors' rapidities, which will be handed
        # to FastJet. In general, this should speed things up a bit.
        if(self.input_vecs_cyl is not None):
            self.fastjet_interface.SetInputsWithEtaPhiM(self.input_vecs.flatten(), self.input_vecs_cyl[:,1:].flatten())
        else:
            self.fastjet_interface.SetInputs(self.input_vecs.flatten()) # will be slower

        # Cluster the jets. They're internally pt-sorted and placed in a map (instead of a vector).
        # TODO: Add user_info to pseudojets. Mostly relevant for ghost-associator.
        self.fastjet_interface.ClusterJets()

        # Previously, we carried around fastjet objects in self.jets_dict. Now we'll avoid doing that explicitly
        jet_indices = np.array(self.fastjet_interface.GetJetIndices())
        self.jet_ordering = np.arange(len(jet_indices))
        self.pt_sorting = self.jet_ordering
        self._jetsToVectors()
        return

    def _jetsToVectors(self):
        """
        This function fetches the jet momenta from our FastjetInterface::JetFinder.
        """
        jet_vectors = self.fastjet_interface.GetJetMomenta() # std::map<Int_t, vector<Double_t>>
        jet_vectors_cyl = self.fastjet_interface.GetJetMomentaCylindrical() # std::map<Int_t, vector<Double_t>>

        self.jet_vectors = {x.first:np.array([*x.second]) for x in jet_vectors} # iterating over the map gives the keys, weirdly enough
        self.jet_vectors_cyl = {x.first:np.array([*x.second]) for x in jet_vectors_cyl}

    @profile_method('FastJetFinderBase._ptSort')
    def _ptSort(self, truncate=False):
        """
        Sorts jets by decreasing pT, and truncates to
        take only the first self.n_jets_max jets.

        This is accomplished by modifying self.jet_ordering.
        """

        if(len(self.jet_vectors) < 1): # no jets -> nothing to do (for case of 1 jet, we'll modify some stuff below)
            return
        elif(len(self.jet_vectors) == 1): # 1 jet -> not much to do, just make sure self.jet_ordering reflects this
            self.jet_ordering = [self.jet_ordering[0]]
            return

        jet_pt = np.array([self.jet_vectors_cyl[i][0] for i in self.jet_ordering])
        is_sorted = np.all(jet_pt[:-1] >= jet_pt[1:])
        needs_truncation = (self.n_jets_max is not None and truncate and len(jet_pt) > self.n_jets_max)

        # Early return if already sorted and no truncation needed
        if is_sorted and not needs_truncation:
            # In principle, we could check that self.pt_sorting is not None,
            # but calling _clusterJets() will ensure this is filled.
            # Still need to set pt_sorting if it doesn't exist
            return

        with profile_block('FastJetFinderBase._ptSort - Sort'):
            if (not is_sorted):
                self.pt_sorting = np.argsort(-jet_pt)

            if(needs_truncation):
                self.pt_sorting = self.pt_sorting[:self.n_jets_max]

            if(len(self.pt_sorting) == 1):
                self.jet_ordering = [self.jet_ordering[self.pt_sorting[0]]]
            else:
                self.jet_ordering = list(operator.itemgetter(*self.pt_sorting)(self.jet_ordering))

            # If we've truncated, we'll need to update some dictionaries under-the-hood to reflect this.
            if(needs_truncation):

                # We will recompute the jet vectors, to account for any that have been dropped.
                # Note that due to the dictionary-based approach, we don't have to recompute this if
                # the jet ordering has simply changed.
                self._jetsToVectors()

                # Also refresh constituents. Again, this only needs to be called if jets were dropped,
                # since its a dictionary so a simple reordering of the jets in self.jet_ordering does
                # not necessitate any change.
                self._fetchJetConstituents()
            return

    @profile_method('FastJetFinderBase._fetchJetConstituents')
    def _fetchJetConstituents(self):
        constituent_indices = self.fastjet_interface.GetJetConstituentIndices() # returns std::map<Int_t, vector<Int_t>>
        self.constituent_indices =     {x.first:list(x.second) for x in constituent_indices} # can iterate over std::map like this, to make dictionary
        self.constituent_vectors =     {i:self.input_vecs[x]           for i,x in self.constituent_indices.items()}
        self.constituent_vectors_cyl = {i:self.input_vecs_cyl[x]       for i,x in self.constituent_indices.items()}

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return