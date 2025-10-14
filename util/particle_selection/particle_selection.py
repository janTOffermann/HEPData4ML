import numpy as np
import ROOT as rt

from util.reconstruction.setup import NtupleProducerSetup
from typing import Optional ,List, TYPE_CHECKING
from util.misc.timing import profile_method


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



# class AlgoSelection(BaseSelector):
#     def __init__(self, selection_algo, n, fixed_length=False):
#         super().__init__()
#         self.particle_selection_algo = selection_algo
#         self.SetN(n)
#         self.fixed_length = fixed_length
#         self.selection_status = True

#     @profile_method('AlgoSelection.__call__')
#     def __call__(self, hepev: 'hm.GenEvent') -> np.ndarray:
#         self.selection_status, particle_list = self.particle_selection_algo(hepev)

#         if self.n > 0 and len(particle_list) > self.n:
#             particle_list = particle_list[:self.n]

#         return particle_list

# class MultiSelection(BaseSelector):
#     def __init__(self, particle_selection_list: List[BaseSelector], enforce_unique=False):
#         super().__init__()
#         self.particle_selection_list = particle_selection_list
#         self.n = sum(x.GetN() for x in particle_selection_list)

#         self.fixed_length = all(x.IsFixedLength() for x in particle_selection_list)
#         self.enforce_unique = enforce_unique
#         self.selection_status = True

#     @profile_method('MultiSelection.__call__')
#     def __call__(self, hepev: 'hm.GenEvent') -> Optional[np.ndarray]:
#         self.selection_status = True
#         particle_lists = []

#         for selector in self.particle_selection_list:
#             particle_list = selector(hepev)

#             if particle_list is None or not selector.GetSelectionStatus():
#                 self.selection_status = False
#                 break

#             # Ensure it's a numpy array
#             if not isinstance(particle_list, np.ndarray):
#                 if np.isscalar(particle_list):
#                     particle_list = np.array([particle_list], dtype=np.int32)
#                 else:
#                     particle_list = np.array(particle_list, dtype=np.int32)

#             particle_lists.append(particle_list)

#         if not self.selection_status or not particle_lists:
#             return None

#         # Concatenate all particle lists
#         result = np.concatenate(particle_lists)

#         if self.enforce_unique:
#             result = np.unique(result)

#         return result