import ROOT as rt
from util.reconstruction.setup import NtupleProducerSetup

class BaseSelectorAlgorithm:
    def __init__(self):
        self.setup = NtupleProducerSetup()
        self.setup.FullPreparation()
        self.algorithm = None

    def GetAlgorithm(self):
        return self.algorithm

class SelectDaughters(BaseSelectorAlgorithm):
    def __init__(self, truth_selection):
        super().__init__()
        self.truth_selection = truth_selection
        self.algorithm = rt.NtupleProducer.SelectDaughters(self.truth_selection.GetSelector())

class SelectFinalStateDaughters(BaseSelectorAlgorithm):
    def __init__(self, truth_selection):
        super().__init__()
        self.truth_selection = truth_selection
        self.algorithm = rt.NtupleProducer.SelectStableDaughters(self.truth_selection.GetSelector())
