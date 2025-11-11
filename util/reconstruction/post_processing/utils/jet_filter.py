import numpy as np
from typing import TYPE_CHECKING
from util.reconstruction.post_processing.utils.postprocessor_base import PostProcessorBase

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder


class PtFilter(PostProcessorBase):
    """
    Removes jets that fail a minimum pT cut (in GeV).
    """

    def __init__(self,pt_min):
        super().__init__()

        self.pt_min = pt_min
        self.name = 'PtFilter'
        self.print_prefix = '{}'.format(self.name)
        self.citations = {}

    def ModifyInitialization(self,obj):
        if(self.obj_name_input is None):
            self.obj_name_input = obj.jet_name
        self.obj_name_output = self.obj_name_input
        return

    def ModifyJets(self, obj:'JetFinder'):
        """
        This function applies the pt filter.
        """

        pmu_cyl_key = '{}.Pmu_cyl'.format(self.obj_name_input)
        pmu_cyl_dict = obj.output_buffer.get(pmu_cyl_key,filter=obj.jet_ordering)

        tags = {i:False for i in pmu_cyl_dict.keys()}

        for i,pmu in pmu_cyl_dict.items():
            if(pmu[0] > self.pt_min):
                tags[i] = True

        obj.jet_ordering = [key for key in obj.jet_ordering if tags[key]]
        # obj._updateJetDictionary()
        # obj._ptSort()
        # obj._jetsToVectors()
        # obj._fetchJetConstituents()
        return

class EtaFilter(PostProcessorBase):
    """
    Removes jets that fail an eta window cut.
    Note: This is a cut on the eta of the jet vector --
          in its current implementation, there isn't a distinction
          between physical eta and detector eta.
    """

    def __init__(self,eta_max):
        super().__init__()

        self.eta_max = eta_max
        self.name = 'EtaFilter'
        self.print_prefix = '{}'.format(self.name)
        self.citations = {}

    def ModifyInitialization(self,obj:'JetFinder'):
        if(self.obj_name_input is None):
            self.obj_name_input = obj.jet_name
        self.obj_name_output = self.obj_name_input

        # Set maximum rapidity for ghosts for jet area calculation;
        # only relevant if using jet areas.
        obj.SetAreaGhostMaxRapidity(self.eta_max + 1.05 * obj.radius)

        return

    def ModifyJets(self, obj:'JetFinder'):
        """
        This function applies the eta filter.
        """
        pmu_cyl_key = '{}.Pmu_cyl'.format(self.obj_name_input)
        pmu_cyl_dict = obj.output_buffer.get(pmu_cyl_key,filter=obj.jet_ordering)

        tags = {i:False for i in pmu_cyl_dict.keys()}

        for i,pmu in pmu_cyl_dict.items():
            if(np.abs(pmu[1]) < self.eta_max):
                tags[i] = True

        obj.jet_ordering = [key for key in obj.jet_ordering if tags[key]]
        # obj._updateJetDictionary()
        # obj._ptSort()
        # obj._jetsToVectors()
        # obj._fetchJetConstituents()
        return

class Leading(PostProcessorBase):
    """
    Removes all jets except the highest-pT one.
    """

    def __init__(self):
        super().__init__()

        self.name = 'Leading'
        self.print_prefix = '{}'.format(self.name)
        self.citations = {}

    def ModifyInitialization(self,obj:'JetFinder'):
        if(self.obj_name_input is None):
            self.obj_name_input = obj.jet_name
        self.obj_name_output = self.obj_name_input

        obj.single_jet = True
        return

    def ModifyInputs(self,obj:'JetFinder'):
        return

    def ModifyJets(self, obj:'JetFinder'):
        """
        This function applies the leading (highest-pT) cut.
        """
        pmu_cyl_key = '{}.Pmu_cyl'.format(self.obj_name_input)
        pmu_cyl_dict = obj.output_buffer.get(pmu_cyl_key,filter=obj.jet_ordering)

        tags = {i:False for i in pmu_cyl_dict.keys()}

        jet_pt = np.array([pmu[0] for pmu in pmu_cyl_dict.values()])
        tags[list(pmu_cyl_dict.keys())[np.argmax(jet_pt)]] = True

        obj.jet_ordering = [key for key in obj.jet_ordering if tags[key]]
        # obj._updateJetDictionary()
        # # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
        # obj._ptSort()
        # obj._jetsToVectors()
        # obj._fetchJetConstituents()
        return