# Some very generic/simple algorithms, which are
# used in selection_algos.py. Keeping them here
# to limit code clutter elsewhere.

import numpy as np
from typing import List, TYPE_CHECKING
from collections import deque

if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports
    from util.hepmc.setup import HepMCSetup, uncache_hepmc3, prepend_to_pythonpath

    # make sure Python HepMC3 bindings are setup
    setup = HepMCSetup(verbose=False)
    # setup.PrepHepMC() # will download/install if necessary
    python_dir = setup.GetPythonDirectory()

    uncache_hepmc3()
    prepend_to_pythonpath(python_dir)

    from pyHepMC3 import HepMC3 as hm

#==========================================
# Here are a bunch of convenience
# functions, which may be useful
# in defining particle selection algorithms.
# By default they take a PythiaWrapper and a
# particle index, but they can also be passed
# a tuple (momentum, status, pdgid) as given
# by PythiaWrapper.GetParticle(), or a
# HepMC::GenParticle (from pyhepmc).
# The priority of inputs is:
# GenParticle > tuple > PythiaWrapper+idx.
#==========================================

def GetDaughtersSingle(hepev: 'hm.GenEvent', idx: int) -> List[int]:
    """Keep original function but with type hints"""
    return [x.id() - 1 for x in hepev.particles()[idx].children()]

def IsStable(hepev: 'hm.GenEvent', idx: int) -> bool:
    """Keep original function"""
    return hepev.particles()[idx].status() == 1

def IsQuark(hepev: 'hm.GenEvent', idx: int) -> bool:
    """Keep original function"""
    pid = abs(hepev.particles()[idx].pid())
    return 0 < pid < 7

class GatherStableDaughters:
    """
    Minimal optimization: Use iterative approach instead of recursion,
    but keep the same basic algorithm structure
    """
    
    def __call__(self, hepev: 'hm.GenEvent', idx: int) -> np.ndarray:
        """Use iterative BFS instead of recursion to avoid function call overhead"""
        stable_daughters = []
        to_visit = deque([idx])
        visited = set([idx])  # Don't include the starting particle
        
        while to_visit:
            current_idx = to_visit.popleft()
            
            # Get daughters using the original method (it might be optimized in HepMC)
            daughters = GetDaughtersSingle(hepev, current_idx)
            
            for daughter_idx in daughters:
                if daughter_idx in visited:
                    continue
                visited.add(daughter_idx)
                
                if IsStable(hepev, daughter_idx):
                    stable_daughters.append(daughter_idx)
                else:
                    to_visit.append(daughter_idx)
        
        return np.unique(np.array(stable_daughters, dtype=np.int32))

class GatherQuarks:
    """Keep closer to original structure but use iterative approach"""
    
    def __call__(self, hepev: 'hm.GenEvent', idx: int) -> np.ndarray:
        quark_daughters = []
        to_visit = deque([idx])
        visited = set([idx])
        
        while to_visit:
            current_idx = to_visit.popleft()
            daughters = GetDaughtersSingle(hepev, current_idx)
            
            for daughter_idx in daughters:
                if daughter_idx in visited:
                    continue
                visited.add(daughter_idx)
                
                if IsQuark(hepev, daughter_idx):
                    if daughter_idx not in quark_daughters:
                        quark_daughters.append(daughter_idx)
                else:
                    to_visit.append(daughter_idx)
        
        return np.array(quark_daughters, dtype=np.int32)

class GatherDaughters:
    """Keep original approach but cleaner"""
    def __init__(self, recursive=False):
        self.recursive = recursive

    def __call__(self, hepev: 'hm.GenEvent', idx: int) -> np.ndarray:
        daughters = GetDaughtersSingle(hepev, idx)
        
        # Handle self-reference (particle re-listing)
        current_particle = hepev.particles()[idx]
        current_pid = current_particle.pid()
        
        for daughter_idx in daughters:
            daughter_particle = hepev.particles()[daughter_idx]
            if daughter_particle.pid() == current_pid:
                # Found self-reference, recurse from there
                return self(hepev, daughter_idx)
        
        return np.array(daughters, dtype=np.int32)