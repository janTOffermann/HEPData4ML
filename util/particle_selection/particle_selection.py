import numpy as np

from typing import Optional ,List, TYPE_CHECKING
from util.misc.timing import profile_method


if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports
    from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath

    # make sure Python HepMC3 bindings are setup
    setup = HepMCSetup(verbose=False)
    python_dir = setup.GetPythonDirectory()

    prepend_to_pythonpath(python_dir)

    from pyHepMC3 import HepMC3 as hm

# ==============================
# These are particle selectors, that the user should
# use (via the configuration). Some are built to
# use specialized algorithms, which are contained
# in selection_algos.py.
# ==============================

class BaseSelector:
    def __init__(self):
        self.selection_status = True
        self.fixed_length = True
        self.n = 1

    def GetN(self):
        return self.n

    def SetN(self, n: int):
        self.n = n

    def GetSelectionStatus(self):
        return self.selection_status

    def IsFixedLength(self):
        return self.fixed_length

class FirstSelector(BaseSelector):
    def __init__(self, status: int, pdgid: int, hadronization: bool = True):
        super().__init__()
        self.status = status
        self.pdgid = pdgid
        self.hadronization = hadronization
        if not hadronization and abs(pdgid) in [1, 2, 3, 4, 5]:
            self.status = 1

    def SetStatus(self, status: int):
        self.status = status

    def SetPdgId(self, pdgid: int):
        self.pdgid = pdgid

    def SetHadronization(self, hadronization: bool = True):
        self.hadronization = hadronization
        if not hadronization and abs(self.pdgid) in [1, 2, 3, 4, 5]:
            self.status = 1

    @profile_method('FirstSelector.__call__')
    def __call__(self, hepev: 'hm.GenEvent') -> Optional[int]:
        particles = hepev.particles()

        pdgid_array = np.array([p.pid() for p in particles], dtype=np.int32)
        status_array = np.array([p.status() for p in particles], dtype=np.int32)

        # Use boolean indexing
        pdgid_mask = pdgid_array == self.pdgid

        if self.status is not None:
            status_mask = status_array == self.status
            combined_mask = pdgid_mask & status_mask
        else:
            combined_mask = pdgid_mask

        indices = np.where(combined_mask)[0]

        if len(indices) == 0:
            self.selection_status = False
            return None

        self.selection_status = True
        return int(indices[0])

    def Print(self):
        print('FirstSelector: status = {}, pdgid = {}'.format(self.status, self.pdgid))
        return

class BasicSelection:
    def __init__(self, selection_list: List[BaseSelector], hadronization=True):
        self.hadronization = hadronization
        self.selection_list = selection_list
        self.n = len(self.selection_list)
        self.fixed_length = True
        self.selection_status = True

    def SetHadronization(self, hadronization=True):
        self.hadronization = hadronization

    @profile_method('BasicSelector.__call__')
    def __call__(self, hepev: 'hm.GenEvent') -> Optional[np.ndarray]:
        self.selection_status = True
        particle_list = []

        for selector in self.selection_list:
            selector.SetHadronization(self.hadronization)
            particle = selector(hepev)
            if particle is not None:
                particle_list.append(particle)
            else:
                self.selection_status = False

        if not particle_list:
            return None

        return np.sort(np.array(particle_list, dtype=np.int32))

class AlgoSelection(BaseSelector):
    def __init__(self, selection_algo, n, fixed_length=False):
        super().__init__()
        self.particle_selection_algo = selection_algo
        self.SetN(n)
        self.fixed_length = fixed_length
        self.selection_status = True

    @profile_method('AlgoSelection.__call__')
    def __call__(self, hepev: 'hm.GenEvent') -> np.ndarray:
        self.selection_status, particle_list = self.particle_selection_algo(hepev)

        if self.n > 0 and len(particle_list) > self.n:
            particle_list = particle_list[:self.n]

        return particle_list

class MultiSelection(BaseSelector):
    def __init__(self, particle_selection_list: List[BaseSelector], enforce_unique=False):
        super().__init__()
        self.particle_selection_list = particle_selection_list
        self.n = sum(x.GetN() for x in particle_selection_list)

        self.fixed_length = all(x.IsFixedLength() for x in particle_selection_list)
        self.enforce_unique = enforce_unique
        self.selection_status = True

    @profile_method('MultiSelection.__call__')
    def __call__(self, hepev: 'hm.GenEvent') -> Optional[np.ndarray]:
        self.selection_status = True
        particle_lists = []

        for selector in self.particle_selection_list:
            particle_list = selector(hepev)

            if particle_list is None or not selector.GetSelectionStatus():
                self.selection_status = False
                break

            # Ensure it's a numpy array
            if not isinstance(particle_list, np.ndarray):
                if np.isscalar(particle_list):
                    particle_list = np.array([particle_list], dtype=np.int32)
                else:
                    particle_list = np.array(particle_list, dtype=np.int32)

            particle_lists.append(particle_list)

        if not self.selection_status or not particle_lists:
            return None

        # Concatenate all particle lists
        result = np.concatenate(particle_lists)

        if self.enforce_unique:
            result = np.unique(result)

        return result