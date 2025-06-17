import sys, operator
import numpy as np
from util.fastjet.setup import FastJetSetup

class JetFinderBase:
    """
    This is a base class, for using the Fastjet library to perform jet clustering.
    It is leveraged by the JetFinder class in utils/post_processing/jets.py .
    This base class exists to avoid some possible import loops.
    It lacks some useful functionality that JetFinder includes.
    """

    def __init__(self,fastjet_dir=None):

        self.jet_algorithm_name = ''
        self.jet_name = ''
        self.radius = 0.4

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
        self.configurator = None

        self.jet_ordering = None

    def _initialize_fastjet(self):

        if(self.fastjet_init_flag):
            return

        if(self.fastjet_dir is None):
            if(self.configurator is None):
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

    def _clusterJets(self):
        self._initialize_fastjet()
        import fastjet as fj # hacky, but will work because _setupFastJet() was run in __init__()

        # vecs has format (E,px,py,pz) -- FastJet uses (px,py,pz,E) so we must modify it. Using np.roll.
        pj = [fj.PseudoJet(*x) for x in np.roll(self.input_vecs,-1,axis=-1)]

        # Attach indices to the pseudojet objects, so that we can trace them through jet clustering.
        # Indices will correspond to the order they were input (with zero-indexing).
        for i,pseudojet in enumerate(pj):
            pseudojet.set_user_index(i)

        # Attach any optional information to the pseudojet objects. This can be leveraged by other classes
        # or extensions.
        if(self.user_info is not None):
            for idx, val in self.user_info.items(): # in practice, val will be a dictionary itself -- allows for attaching multiple things
                pj[idx].set_python_info(val)


        # selector = fj.SelectorPtMin(jet_config['jet_min_pt']) & fj.SelectorAbsEtaMax(jet_config['jet_max_eta'])
        # Note: Switched from the old method, this is more verbose but seems to do the same thing anyway.
        self.cluster_sequence = fj.ClusterSequence(pj, self.jetdef) # member of class, otherwise goes out-of-scope when ref'd later
        # self.jets = self.cluster_sequence.inclusive_jets()
        self.jets_dict = {i:jet for i,jet in enumerate(self.cluster_sequence.inclusive_jets())}
        self.jet_ordering = np.arange(len(self.jets_dict))
        self._jetsToVectors()

    def _jetsToVectors(self):
        """
        This function fills self.jet_vectors and self.jet_vectors_cyl, to contain the four-momenta
        of whatever jets are currently in self.jets.
        """
        # self.jet_vectors = np.array([[x.e(), x.px(), x.py(),x.pz()] for x in self.jets])
        # self.jet_vectors_cyl = np.array([[x.pt(), x.eta(), x.phi(),x.m()] for x in self.jets])

        self.jet_vectors = {i:np.array([jet.e(), jet.px(), jet.py(),jet.pz()]) for i,jet in self.jets_dict.items()}
        self.jet_vectors_cyl = {i:np.array([jet.pt(), jet.eta(), jet.phi(),jet.m()]) for i,jet in self.jets_dict.items()}

    def _ptSort(self, truncate=False):
        """
        Sorts jets by decreasing pT, and truncates to
        take only the first self.n_jets_max jets.

        This is accomplished by modifying self.jet_ordering.
        """

        if(len(self.jets_dict) < 1): # no jets -> nothing to do
            return

        jet_pt = np.array([self.jets_dict[i].pt() for i in self.jet_ordering])
        self.pt_sorting = np.argsort(-jet_pt)

        if(self.n_jets_max is not None):
            if((len(self.pt_sorting) > self.n_jets_max) and truncate):
                self.pt_sorting = self.pt_sorting[:self.n_jets_max]

        if(len(self.pt_sorting) == 1):
            self.jet_ordering = [self.jet_ordering[self.pt_sorting[0]]]
        else:
            self.jet_ordering = list(operator.itemgetter(*self.pt_sorting)(self.jet_ordering))

        if(len(self.jet_ordering) != len(self.jets_dict)):
            # remove entries from jets_dict, that correspond with entries in jet_ordering that have been dropped
            self._updateJetDictionary()

            # We will also recompute the jet vectors, to account for any that have been dropped.
            # Note that due to the dictionary-based approach, we don't have to recompute this if
            # the jet ordering has simply changed.
            self._jetsToVectors()

            # Also refresh constituents.
            self._fetchJetConstituents()

        self._fetchJetConstituents()

    def _updateJetDictionary(self):
        """
        To be used when self.jet_ordering is updated.
        """
        for key in list(self.jets_dict.keys()):
            if(key not in self.jet_ordering):
                del self.jets_dict[key]
        return

    def _fetchJetConstituents(self):
        results = {i:self._fetchJetConstituentsSingle(jet, self.n_constituents_max) for i,jet in self.jets_dict.items()}
        self.constituent_vectors = {i:x[0] for i,x in results.items()}
        self.constituent_vectors_cyl = {i:x[1] for i,x in results.items()}
        self.constituent_indices = {i:x[2] for i,x in results.items()}

    def _fetchJetConstituentsSingle(self,jet,n_constituents=-1):

        # TODO: Understand if/when this is called on a jet without constituents (i.e. a base PseudoJet).
        if(not jet.has_constituents()):
            return np.empty((0, 4)), np.empty((0, 4)), np.empty((0, 4)) # TODO: Is this OK?

        pt,eta,phi,m,e,px,py,pz = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), x.m(), x.e(), x.px(),x.py(),x.pz()] for x in jet.constituents()]),8)

        # The indices of the jet constituents, corresponding with the order in which they
        # were passed to jet clustering.
        indices = np.array([x.user_index() for x in jet.constituents()],dtype=np.dtype('i4')).flatten()

        # Sort by decreasing pt, and only keep leading constituents.
        sorting = np.argsort(-pt.flatten())
        l = len(pt)
        if((n_constituents is not None) and n_constituents > 0):
            l = int(np.minimum(n_constituents,l))
        vecs = np.vstack([x.flatten() for x in [e,px,py,pz]]).T[sorting][:l]
        vecs_cyl = np.vstack([x.flatten() for x in [pt,eta,phi,m]]).T[sorting][:l]
        indices = indices[sorting][:l]

        return vecs, vecs_cyl, indices

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