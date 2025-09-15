import numpy as np
from util.particle_selection.algos import *
from typing import Union, TYPE_CHECKING
from util.misc.timing import profile_method, profile_block

if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports
    from util.hepmc.setup import HepMCSetup, uncache_hepmc3, prepend_to_pythonpath

    # make sure Python HepMC3 bindings are setup
    setup = HepMCSetup(verbose=False)
    # setup.PrepHepMC() # will download/install if necessary
    python_dir = setup.GetPythonDirectory()

    # uncache_hepmc3()
    prepend_to_pythonpath(python_dir)

    from pyHepMC3 import HepMC3 as hm

# ============================================
# Each algorithm returns a Boolean status,
# and a list of selected particles' event listing indices.
# These can be used by AlgoSelection() in particle_selection.py.
# ============================================

class BaseSelectorAlgorithm:
    """Keep simple"""
    def __init__(self):
        pass

class SelectFinalState(BaseSelectorAlgorithm):
    """Minimal optimization: vectorize the status check"""
    def __init__(self):
        pass

    @profile_method('SelectFinalState.__call__')
    def __call__(self, hepev: 'hm.GenEvent') -> tuple:
        particles = hepev.particles()
        status_array = np.array([p.status() for p in particles], dtype=np.int32)
        stable_indices = np.where(status_array == 1)[0]
        return len(stable_indices) > 0, stable_indices

class SelectDaughters(BaseSelectorAlgorithm):
    """Keep original structure but use optimized gatherer"""
    def __init__(self, truth_selection):
        self.truth_selection = truth_selection
        self.gatherer = GatherDaughters()

    @profile_method('SelectDaughters.__call__')
    def __call__(self, hepev: 'hm.GenEvent') -> tuple:
        starting_particles = self.truth_selection(hepev)
        if not isinstance(starting_particles, (list, np.ndarray)):
            starting_particles = [starting_particles]
        
        all_daughters = []
        for p in starting_particles:
            daughters = self.gatherer(hepev, p)
            all_daughters.extend(daughters)
        
        if not all_daughters:
            return False, np.array([], dtype=np.int32)

        return True, np.unique(np.array(all_daughters, dtype=np.int32))

class SelectFinalStateDaughters(BaseSelectorAlgorithm):
    """Use the optimized gatherer but keep simple interface"""
    def __init__(self, truth_selection):
        self.truth_selection = truth_selection
        self.gatherer = GatherStableDaughters()

    @profile_method('SelectFinalStateDaughters.__call__')
    def __call__(self, hepev: 'hm.GenEvent') -> tuple:
        starting_particles = self.truth_selection(hepev)
        if not isinstance(starting_particles, (list, np.ndarray)):
            starting_particles = [starting_particles]
        
        if not starting_particles or starting_particles[0] is None:
            return False, np.array([], dtype=np.int32)
        
        all_daughters = []
        for start_idx in starting_particles:
            if start_idx is None or start_idx < 0:
                continue
            daughters = self.gatherer(hepev, start_idx)
            all_daughters.extend(daughters)
        
        if not all_daughters:
            return False, np.array([], dtype=np.int32)
        
        # Use numpy unique for deduplication and sorting
        result = np.unique(np.array(all_daughters, dtype=np.int32))
        return True, result