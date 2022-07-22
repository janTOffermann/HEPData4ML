import numpythia as npyth # for defining selections
import numpy as np
from util.hepmc import HepMCArrayToNumpy, RestructureParticleArray
# from util.truth_selection import selections as truth_sel

def IsQuark(particle):
    pid = np.abs(particle.pid)
    return (pid > 0 and pid < 7)

def IsLepton(particle):
    pid = np.abs(particle.pid)
    return (pid > 10 and pid < 19)

def IsGluon(particle):
    pid = np.abs(particle.pid)
    return (pid in [9,21])

def GatherQuarks(particle,ignore=False):
    particle_list = []
    if(IsQuark(particle) and not ignore):
        particle_list.append(particle)
    elif(not IsGluon(particle)):
        children = particle.children(return_hepmc=True)
        if(len(children) != 0):
            for child in children:
                particle_list += GatherQuarks(child)
    return particle_list

# -----------------------------------------------------------------

# This function just selects all actual final-state (i.e. stable) particles.
def SelectFinalState(**kwargs):
    event = kwargs['event']
    selection = (npyth.STATUS == 1) & (npyth.ABS_PDG_ID > -1)
    particles = event.all(selection=selection, return_hepmc=False)
    return particles

# Select the "simplest" all-quark final state -- this recursively searches the Feynman diagram starting with
# the particles given by "truth_sel", and selects the earliest quark found along each branch.
# This is a callable class instead of a function, so that we can configure it with its own selection!
class SelectSimplestHadronic:
    def __init__(self,truth_sel):
        self.truth_sel = truth_sel

    def __call__(self, event):
        truth_arr = [event.first(selection=x, return_hepmc=True) for x in self.truth_sel]
        particles = []
        for par in truth_arr:
            particles += GatherQuarks(par)
        return HepMCArrayToNumpy(particles)