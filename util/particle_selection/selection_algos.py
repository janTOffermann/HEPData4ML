import numpy as np
from util.particle_selection.algos import *

# ============================================
# Each algorithm returns a Boolean status,
# and a list of selected particles' event listing indices.
# These can be used by AlgoSelection() in particle_selection.py.
# ============================================

# This function just selects all actual final-state (i.e. stable) particles in an event.
class SelectFinalState:
    def __init__(self):
        pass

    def __call__(self,pythia_wrapper):
        status = pythia_wrapper.GetStatus(hepmc=True)
        stable = np.where(status==1)[0]
        return (len(stable) > 0), stable

# Select the "simplest" all-quark final state produced by the particles in truth_selection.
# This recursively searches the Feynman diagram starting with the particles
# given by "truth_sel", and selects the earliest quark found along each branch.
class SelectSimplestQuarks:
    def __init__(self,truth_selection):
        self.truth_selection = truth_selection

    def __call__(self,pythia_wrapper):
        starting_particles = self.truth_selection(pythia_wrapper)
        if(type(starting_particles) != list): starting_particles = [starting_particles]
        particles = []
        for p in starting_particles:
            gatherer = GatherQuarks()
            particles.append(gatherer(pythia_wrapper,p))
        if(len(particles) == 0): return False, particles

        particles = np.array(particles,dtype=int).flatten()
        return True, np.unique(particles)

# Select the stable daughters of the particles chosen by truth_selection.
class SelectFinalStateDaughters:
    def __init__(self,truth_selection):
        self.truth_selection = truth_selection

    def __call__(self,pythia_wrapper):
        starting_particles = self.truth_selection(pythia_wrapper)
        if(type(starting_particles) not in [list,np.ndarray]): starting_particles = [starting_particles]
        particles = []
        for p in starting_particles:
            gatherer = GatherStableDaughters()
            daughters = gatherer(pythia_wrapper,p)
            particles.append(daughters)
        if(len(particles) == 0): return False, particles

        particles = np.concatenate(particles,dtype=int).flatten()
        return True, np.unique(particles)
