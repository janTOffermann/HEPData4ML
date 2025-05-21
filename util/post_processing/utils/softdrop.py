import sys
import numpy as np
from util.fastjet import JetFinderBase
from util.calcs import embed_array_inplace

class Softdrop:
    """
    Performs the Softdrop algorithm.
    See: https://arxiv.org/abs/1402.2657 [JHEP 05 (2014) 146]

    Note: Keep in mind that based on how JetFinder is structured,
          the soft-dropped jets will be pT-sorted *after* the softdrop.
          Thus if you have a separate JetFinder producing the same jets
          but without softdrop, they might not line up jet-by-jet, because
          performing softdrop will change the jet pT and thus possibly ordering.
    """

    def __init__(self,zcut,beta):

        self.zcut = zcut
        self.beta = beta
        self.jet_finder = JetFinderBase()

        self.cluster_sequences = []

    def _initialize_fastjet(self): #TODO: This can probably be removed, as in practice JetFinder will have taken care of this
        """
        Based on JetFinderBase._initialize_fastjet().
        Not sure if this is really needed in practice.
        """
        self.jet_finder._setupFastJet()
        sys.path.append(self.jet_finder.fastjet_dir)
        import fastjet as fj

    def ModifyInitialization(self,obj):
        # Take the opportunity to set up FastJet.
        # In practice, this shouldn't be needed
        # based on how this is called, I think?

        self.jet_finder.SetConfigurator(obj.configurator)
        self.jet_finder.SetRadius(obj.radius)
        self.jet_finder.jet_algorithm_name = 'c/a'

        # Initialize fastjet.
        self.jet_finder._initialize_fastjet()

        self._initialize_fastjet() # TODO: Not sure if this is needed?

        # Initialize the fastjet jet definition
        self.jet_finder._initialize_jet_definition()

        return

    def ModifyInputs(self,obj):
        return

    def ModifyJets(self, obj):
        """
        This function will perform softdrop.
        This involves re-clustering jets with the C/A algorithm.
        """
        # self.old_jets = [] # prevent the pre-softdrop jets from immediately going out-of-scope
        new_jets = []

        # We must store cluster sequences, otherwise they will go out-of-scope.
        # This is an issue if we try to access jet constituents later.
        self.cluster_sequences = []

        for jet in obj.jets:

            # Recluster with C/A
            ca_jet = self._recluster(jet)

            # Now run softdrop
            groomed_jet = self._softdrop(ca_jet)
            new_jets.append(groomed_jet)

        obj.jets = new_jets
        obj._jetsToVectors()

        return

    def ModifyWrite(self,obj):
        return # does nothing

    def _recluster(self,jet):
        # Gather the jet constituents, recluster with C/A
        constituents, _, _ = self.jet_finder._fetchJetConstituentsSingle(jet)
        self.jet_finder.SetInputs(constituents)
        self.jet_finder._clusterJets()

        # Store the cluster sequence from jet finder.
        # Need these to not go out-of-scope if constituents are accessed later.
        # This will get refreshed next time self.ModifyJets() is called.
        # (Python memory management is sometimes quite a pain compared to C++!)
        self.cluster_sequences.append(self.jet_finder.cluster_sequence)

        #NOTE: Based on setup, there should only be one jet. Nonetheless, as a precaution
        #      we will take the leading (highest-pT) one.
        jet_pt = [x.pt for x in self.jet_finder.jets]
        ca_jet = self.jet_finder.jets[np.argmax(jet_pt)]
        return ca_jet

    def _softdrop(self,jet):
        import fastjet as fj

        R = self.jet_finder.radius

        current_jet = jet
        groomed_jet = None
        parent1 = fj.PseudoJet()
        parent2 = fj.PseudoJet()

        while(current_jet.has_parents(parent1,parent2)):

            if(parent2.pt() > parent1.pt()):
                parent1, parent2 = parent2, parent1

            z = parent2.pt() / (parent1.pt() + parent2.pt())
            dR = parent1.delta_R(parent2) # NOTE: This function, from fastjet, computes dR using rapidity, not pseudo-rapidity!
            threshold = self.zcut * np.power(dR/R, self.beta)

            if(z > threshold):
                groomed_jet = current_jet
                break
            else:
                current_jet = parent1 # run again, on the higher-pT jet

        if(groomed_jet is None):
            groomed_jet = current_jet # will be whatever the last hard subjet was

        return groomed_jet

    def ModifyConstituents(self, obj):
        return


class IteratedSoftdrop:
    """
    Performs the Iterated Softdrop (ISD) algorithm.
    See: https://arxiv.org/abs/1704.06266 [JHEP 09 (2017) 083]
    """

    def __init__(self,zcut,beta, dR_cut, max_depth=5, mode='prune',name='IteratedSoftDrop'):

        self.zcut = zcut
        self.beta = beta
        self.dR_cut = dR_cut
        self.jet_finder = JetFinderBase()

        self.max_depth = max_depth
        self.mode = mode.lower()
        assert self.mode in ['prune','tag']

        self.name = name
        self.branch_name_prefix = None
        self.zi_name = None
        self.dRi_name = None

        self.zi = None
        self.dRi = None

        self.cluster_sequences = []

    def ModifyInitialization(self,obj):
        # Take the opportunity to set up FastJet.
        # In practice, this shouldn't be needed
        # based on how this is called, I think?

        self.jet_finder.SetConfigurator(obj.configurator)
        self.jet_finder.SetRadius(obj.radius)
        self.jet_finder.jet_algorithm_name = 'c/a'

        # Initialize fastjet.
        self.jet_finder._initialize_fastjet()

        # Initialize the fastjet jet definition
        self.jet_finder._initialize_jet_definition()

        return

    def ModifyInputs(self,obj):
        return

    def ModifyJets(self, obj):
        """
        This function will perform softdrop.
        This involves re-clustering jets with the C/A algorithm.

        Depending on the run mode, this may actually modify the jets,
        or it may simply compute the ISD observables.
        """
        new_jets = []
        self.zi =  np.zeros((obj.n_jets_max,self.max_depth))
        self.dRi = np.zeros((obj.n_jets_max,self.max_depth))

        # We must store cluster sequences, otherwise they will go out-of-scope.
        # This is an issue if we try to access jet constituents later.
        self.cluster_sequences = []

        for i,jet in enumerate(obj.jets):

            # Recluster with C/A
            ca_jet = self._recluster(jet)

            # Now run softdrop
            groomed_jet = self._softdrop(ca_jet,i) # also fills self.zi, self.dRi for this jet
            new_jets.append(groomed_jet)

        if(self.mode=='prune'):
            obj.jets = new_jets
            obj._jetsToVectors()

        return

    def ModifyWrite(self,obj):
        # always runs
        self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
        self._addToBuffer(obj)

    def _recluster(self,jet):
        # Gather the jet constituents, recluster with C/A
        constituents, _, _ = self.jet_finder._fetchJetConstituentsSingle(jet)
        self.jet_finder.SetInputs(constituents)
        self.jet_finder._clusterJets()

        # Store the cluster sequence from jet finder.
        # Need these to not go out-of-scope if constituents are accessed later.
        # This will get refreshed next time self.ModifyJets() is called.
        # (Python memory management is sometimes quite a pain compared to C++!)
        self.cluster_sequences.append(self.jet_finder.cluster_sequence)

        #NOTE: Based on setup, there should only be one jet. Nonetheless, as a precaution
        #      we will take the leading (highest-pT) one.
        jet_pt = [x.pt for x in self.jet_finder.jets]
        ca_jet = self.jet_finder.jets[np.argmax(jet_pt)]
        return ca_jet

    def _softdrop(self,jet,i):
        import fastjet as fj

        R = self.jet_finder.radius
        N = 0

        current_jet = jet
        groomed_jet = None
        parent1 = fj.PseudoJet()
        parent2 = fj.PseudoJet()

        while(current_jet.has_parents(parent1,parent2)):

            if(parent2.pt() > parent1.pt()):
                parent1, parent2 = parent2, parent1

            z = parent2.pt() / (parent1.pt() + parent2.pt())
            dR = parent1.delta_R(parent2) # NOTE: This function, from fastjet, computes dR using rapidity, not pseudo-rapidity!

            if(dR < self.dR_cut):
                break

            threshold = self.zcut * np.power(dR/R, self.beta)

            if(z > threshold):
                groomed_jet = current_jet

                # save the Zi and dRi
                if(N < self.max_depth):
                    self.zi[i,N] = z
                    self.dRi[i,N] = dR

                # iterate the counter
                N += 1

            current_jet = parent1 # run again, on the higher-pT jet

        if(groomed_jet is None):
            groomed_jet = current_jet # will be whatever the last hard subjet was

        return groomed_jet

    def ModifyConstituents(self, obj):
        return

    def _initializeBuffer(self,obj):
        """
        Adds branches to the buffer corresponding with the Zi and dRi variables
        from the ISD algorithm.
        """
        self._createBranchNames(obj)

        if(self.zi_name not in obj.buffer.keys()):
            obj.buffer[self.zi_name] = np.zeros((obj.nevents,obj.n_jets_max,self.max_depth),dtype=np.dtype('f8'))

        if(self.dRi_name not in obj.buffer.keys()):
            obj.buffer[self.dRi_name] = np.zeros((obj.nevents,obj.n_jets_max,self.max_depth),dtype=np.dtype('f8'))
        return

    def _createBranchNames(self,obj):
        if(self.name is None):
            self.name = 'IteratedSoftdrop'
        self.branch_name_prefixname = '{}.{}'.format(obj.jet_name,self.name)

        self.zi_name = '{}.Z'.format(self.branch_name_prefixname)
        self.dRi_name = '{}.dR'.format(self.branch_name_prefixname)
        return

    def _addToBuffer(self,obj):
        """
        Adds the Zi and dRi info to the buffer.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        #NOTE: The embed is not needed, since we've constructed the inputs and the buffer to already match in size.
        #      The zero-padding is actually being handled within self.ModifyJets(), where the embed function is used.
        embed_array_inplace(self.zi[obj.pt_sorting],obj.buffer[self.zi_name][obj._i])
        embed_array_inplace(self.dRi[obj.pt_sorting],obj.buffer[self.dRi_name][obj._i])

