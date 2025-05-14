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

    def __call__(self,hepev):
        status = [x.status for x in hepev.particles]
        stable = np.where(status==1)[0]
        return (len(stable) > 0), stable

# Select the immediate daughters of some particle.
class SelectDaughters:
    """
    Select the immediate daughter(s) of some particle, in event listing.
    Does not include radiation like X -> X gamma (here, {X,gamma} are *not*
    daughters, we look further at how X decays to find them).
    """
    def __init__(self,truth_selection):
        self.truth_selection = truth_selection
        self.gatherer = GatherDaughters()

    def __call__(self,hepev):
        starting_particles = self.truth_selection(hepev)
        if(type(starting_particles) != list): starting_particles = [starting_particles]
        particles = []
        for p in starting_particles:
            particles.append(self.gatherer(hepev,p))
        if(len(particles) == 0): return False, particles

        particles = np.array(particles,dtype=int).flatten()
        return True, np.unique(particles)

class SelectSimplestQuarks(SelectDaughters):
    """
    Select the "simplest" all-quark final state produced by the particles in truth_selection.
    This recursively searches the Feynman diagram starting with the particles
    given by "truth_sel", and selects the earliest quark found along each branch.
    """
    def __init__(self,truth_selection):
        self.truth_selection = truth_selection
        self.gatherer = GatherQuarks()

# # Select the "simplest" all-quark final state produced by the particles in truth_selection.
# # This recursively searches the Feynman diagram starting with the particles
# # given by "truth_sel", and selects the earliest quark found along each branch.
# class SelectSimplestQuarks:
#     def __init__(self,truth_selection):
#         self.truth_selection = truth_selection

#     def __call__(self,pythia_wrapper):
#         starting_particles = self.truth_selection(pythia_wrapper)
#         if(type(starting_particles) != list): starting_particles = [starting_particles]
#         particles = []
#         for p in starting_particles:
#             gatherer = GatherQuarks()
#             particles.append(gatherer(pythia_wrapper,p))
#         if(len(particles) == 0): return False, particles
#         particles = np.array(particles,dtype=int).flatten()
#         return True, np.unique(particles)

# Select the stable daughters of the particles chosen by truth_selection.
class SelectFinalStateDaughters:
    def __init__(self,truth_selection):
        self.truth_selection = truth_selection

    def __call__(self,hepev):
        starting_particles = self.truth_selection(hepev)
        if(type(starting_particles) not in [list,np.ndarray]): starting_particles = [starting_particles]
        particles = []
        for i,p in enumerate(starting_particles):
            gatherer = GatherStableDaughters()
            daughters = gatherer(hepev,p)
            particles.append(daughters)
        if(len(particles) == 0): return False, particles

        particles = np.unique(np.concatenate(particles,dtype=int).flatten())
        return True, particles
