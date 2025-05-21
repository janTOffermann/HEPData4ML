import sys
import numpy as np
from util.fastjet import JetFinderBase
from util.calcs import embed_array_inplace

class PtFilter:
    """
    Removes jets that fail a minimum pT cut (in GeV).
    """

    def __init__(self,pt_min):
        self.pt_min = pt_min
        self.print_prefix = '\n\t\tPtFilter'

    def ModifyInitialization(self,obj):
        return

    def ModifyInputs(self,obj):
        return

    def ModifyJets(self, obj):
        """
        This function apply the pt filter.
        """
        # self.old_jets = [] # prevent the pre-softdrop jets from immediately going out-of-scope
        new_jets = []

        for jet in obj.jets:
            if(jet.pt() > self.pt_min):
                new_jets.append(jet)

        obj.jets = new_jets
        obj._ptSort()
        # obj._jetsToVectors()

        return

    def ModifyWrite(self,obj):
        return # does nothing

    def ModifyConstituents(self, obj):
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

    def ModifyInitialization(self,obj):
        return

    def ModifyInputs(self,obj):
        return

    def ModifyJets(self, obj):
        """
        This function apply the eta filter.
        """
        # self.old_jets = [] # prevent the pre-softdrop jets from immediately going out-of-scope
        new_jets = []

        for jet in obj.jets:
            if(np.abs(jet.eta()) < self.eta_max):
                new_jets.append(jet)

        obj.jets = new_jets
        obj._ptSort()
        # obj._jetsToVectors()

        return

    def ModifyWrite(self,obj):
        return # does nothing

    def ModifyConstituents(self, obj):
        return

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return
