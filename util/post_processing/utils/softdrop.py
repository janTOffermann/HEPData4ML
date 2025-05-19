import sys
import numpy as np
from util.fastjet import JetFinderBase

class Softdrop():
    """
    Performs the Softdrop algorithm.
    See: https://arxiv.org/abs/1402.2657 [JHEP 05 (2014) 146]

    Note: Keep in mind that based on how JetFinder is structured,
          the soft-dropped jets will be pT-sorted *after* the softdrop.
          Thus if you have another post-processor producing the same jets
          but without softdrop, they might not line up jet-by-jet, because
          performing softdrop will change the jet pT and thus possibly ordering.
    """

    def __init__(self,zcut,beta):

        self.zcut = zcut
        self.beta = beta
        self.jet_finder = JetFinderBase()

        self.cluster_sequences = []

    def _initialize_fastjet(self):
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
            # cluster_sequence = ca_jet.validated_cluster_sequence()
            # self.cluster_sequences.append(cluster_sequence)

            # Now run softdrop
            groomed_jet = self._softdrop(ca_jet)
            new_jets.append(groomed_jet)

        obj.jets = new_jets
        return

    def _recluster(self,jet):
        # Gather the jet constituents, recluster with C/A
        constituents, _, _ = self.jet_finder._fetchJetConstituentsSingle(jet)
        self.jet_finder.SetInputs(constituents)
        self.jet_finder._clusterJets()

        # Store the cluster sequence from jet finder.
        # Need these to not go out-of-scope if constituents are accessed later.
        # This will get refreshed next time self.ModifyJets() is called.
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

            # parent1, parent2 = current_jet.parents(cluster_sequence) # NOTE: assuming length = 2
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

        # # Some hackery, to deal with cluster sequences going out-of-scope.
        # constituents = groomed_jet.constituents()
        # groomed_jet = fj.join(constituents)

        return groomed_jet

    def ModifyConstituents(self, obj):
        return