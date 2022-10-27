import numpythia as npyth # for defining selections
import numpy as np
from util.hepmc import HepMCArrayToNumpy
from util.particle_selection.algos import GatherQuarks, GatherStableDaughters

# This function just selects all actual final-state (i.e. stable) particles.
def SelectFinalState(**kwargs):
    event = kwargs['event']
    selection = (npyth.STATUS == 1) & (npyth.ABS_PDG_ID > -1)
    particles = event.all(selection=selection, return_hepmc=False)
    status = len(particles) != 0
    return status, particles

# Select the "simplest" all-quark final state -- this recursively searches the Feynman diagram starting with
# the particles given by "truth_sel", and selects the earliest quark found along each branch.
# This is a callable class instead of a function, so that we can configure it with its own selection!
class SelectSimplestHadronic:
    def __init__(self,truth_sel, return_hepmc=False):
        self.truth_sel = truth_sel
        self.return_hepmc = return_hepmc

    def SetHepMC(self,flag):
        self.return_hepmc = flag

    def __call__(self, event):
        truth_arr = self.truth_sel(event,return_hepmc=True)
        particles = []
        for par in truth_arr:
            if(par is None): return False, np.zeros(0) # part of truth_sel returned no particle -> this is bad, return status=False (failure)
            particles += GatherQuarks(par)
        if(not self.return_hepmc): particles = HepMCArrayToNumpy(particles)
        return True, particles

class SelectFinalStateDaughters:
    def __init__(self,truth_sel, return_hepmc=False):
        self.truth_sel = truth_sel
        self.return_hepmc = return_hepmc

    def SetHepMC(self,flag):
        self.return_hepmc = flag

    def __call__(self, event): # TODO: Rework this.
        truth_arr = self.truth_sel(event,return_hepmc=True)
        particles = []
        for par in truth_arr:
            if(par is None): return False, np.zeros(0) # part of truth_sel returned no particle -> this is bad, return status=False (failure)
            particles += GatherStableDaughters(par)
        if(not self.return_hepmc): particles = HepMCArrayToNumpy(particles)
        return True, particles
