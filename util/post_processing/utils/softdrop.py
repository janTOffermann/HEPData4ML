
import numpy as np
import itertools
from util.fastjet import JetFinderBase

class Softdrop():

    def __init__(self):

        self.jet_finder = JetFinderBase()
        pass

    def ModifyInitialization(self,obj):
        # Take the opportunity to set up FastJet.
        # In practice, this shouldn't be needed
        # based on how this is called, I think?

        self.jet_finder.SetConfigurator(obj.configurator)
        self.jet_finder.SetRadius(obj.radius)
        self.jet_finder.jet_algorithm_name = 'kt' # TODO: Change to Cambridge/Aachen

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
        """

        new_jets = []

        for jet in obj.jets:
            pass

        # # Because of the indexing done in ModifyInitialization(), here
        # # we know that the ghosts correspond to th first N particles.
        # indices = np.arange(obj.input_collection_arrays[self.ghost_key].shape[1])

        # # Fetch the jet constituents.
        # obj._fetchJetConstituents()

        # mask = np.full(len(obj.constituent_indices),False)

        # for i,constituent_list in enumerate(obj.constituent_indices):
        #     for idx in indices:
        #         if(idx in constituent_list):
        #             mask[i] = True
        #             break

        # obj.jets = list(itertools.compress(obj.jets,mask)) # NOTE: Not really necessary
        # obj.jet_vectors = list(itertools.compress(obj.jet_vectors,mask))
        # obj.jet_vectors_cyl = list(itertools.compress(obj.jet_vectors_cyl,mask))

        return

    def ModifyConstituents(self, obj):
        return