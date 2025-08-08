import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.post_processing.jets import JetFinder


class PtFilter:
    """
    Removes jets that fail a minimum pT cut (in GeV).
    """

    def __init__(self,pt_min):
        self.pt_min = pt_min
        self.print_prefix = '\n\t\tPtFilter'

    def ModifyInitialization(self,obj):
        return

    def ModifyInputs(self,obj:'JetFinder'):
        return

    def ModifyJets(self, obj:'JetFinder'):
        """
        This function apply the pt filter.
        """
        tags = {i:False for i in obj.jets_dict.keys()}

        for i,jet in obj.jets_dict.items():
            if(jet.pt() > self.pt_min):
                tags[i] = True

        obj.jet_ordering = [key for key in obj.jet_ordering if tags[key]]
        obj._updateJetDictionary()
        # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
        obj._ptSort()
        obj._jetsToVectors()
        obj._fetchJetConstituents()
        return

    def ModifyWrite(self,obj:'JetFinder'):
        return # does nothing

    def ModifyConstituents(self, obj:'JetFinder'):
        return

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return

class EtaFilter:
    """
    Removes jets that fail an eta window cut.
    Note: This is a cut on the eta of the jet vector --
          in its current implementation, there isn't a distinction
          between physical eta and detector eta.
    """

    def __init__(self,eta_max):
        self.eta_max = eta_max
        self.print_prefix = '\n\t\tEtaFilter'

    def ModifyInitialization(self,obj:'JetFinder'):
        return

    def ModifyInputs(self,obj):
        return

    def ModifyJets(self, obj):
        """
        This function apply the eta filter.
        """
        tags = {i:False for i in obj.jets_dict.keys()}

        for i,jet in obj.jets_dict.items():
            if(np.abs(jet.eta()) < self.eta_max):
                tags[i] = True

        obj.jet_ordering = [key for key in obj.jet_ordering if tags[key]]
        obj._updateJetDictionary()
        # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
        obj._ptSort()
        obj._jetsToVectors()
        obj._fetchJetConstituents()
        return

    def ModifyWrite(self,obj:'JetFinder'):
        return # does nothing

    def ModifyConstituents(self, obj:'JetFinder'):
        return

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return

class Leading:
    """
    Removes all jets except the highest-pT one.
    """

    def __init__(self):
        self.print_prefix = '\n\t\tLeading'

    def ModifyInitialization(self,obj:'JetFinder'):
        obj.single_jet = True
        return

    def ModifyInputs(self,obj:'JetFinder'):
        return

    def ModifyJets(self, obj:'JetFinder'):
        """
        This function apply the leading (highest-pT) cut.
        """
        tags = {i:False for i in obj.jets_dict.keys()}

        jet_pt = np.array([jet.pt() for jet in obj.jets_dict.values()])
        tags[list(obj.jets_dict.keys())[np.argmax(jet_pt)]] = True

        obj.jet_ordering = [key for key in obj.jet_ordering if tags[key]]
        obj._updateJetDictionary()
        # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
        obj._ptSort()
        obj._jetsToVectors()
        obj._fetchJetConstituents()
        return

    def ModifyWrite(self,obj:'JetFinder'):
        return # does nothing

    def ModifyConstituents(self, obj:'JetFinder'):
        return

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return