import ROOT as rt

from util.reconstruction.setup import NtupleProducerSetup
from typing import List, TYPE_CHECKING


if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports
    from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath

    # make sure Python HepMC3 bindings are setup
    setup = HepMCSetup(verbose=False)
    python_dir = setup.GetPythonDirectory()

    prepend_to_pythonpath(python_dir)

    from pyHepMC3 import HepMC3 as hm

# ==========================================#
# Wrapper classes for particle selectors;   #
# this lets us hide the setup of the        #
# NtupleProducer library, so it's not       #
# sitting in the configuration file itself. #
# ==========================================#


class BaseSelector:
    def __init__(self):

        self.setup = NtupleProducerSetup()
        self.setup.FullPreparation()

        self.selector = None

    def GetSelector(self):
        return self.selector # actual access to the underlying ROOT.NtupleProducer selector

class FirstSelector(BaseSelector):
    def __init__(self, status: int, pdgid: int, hadronization: bool = True):
        super().__init__()
        self.selector = rt.NtupleProducer.FirstSelector(status,pdgid,hadronization)

class AlgoSelection(BaseSelector):
    def __init__(self, algorithm, n=-1):
        super().__init__()
        self.selector = rt.NtupleProducer.AlgoSelection(algorithm.GetAlgorithm(),n)

class MultiSelection(BaseSelector):
    def __init__(self, particle_selection_list:List[BaseSelector], enforce_unique:bool=False):
        super().__init__()
        self.selection_list = [x.GetSelector() for x in particle_selection_list]
        self.selector = rt.NtupleProducer.MultiSelection(self.selection_list, enforce_unique)