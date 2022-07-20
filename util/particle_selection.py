import numpythia as npyth # for defining selections
import numpy as np
from util.truth_selection import selections as truth_sel

def IsQuark(particle):
    pid = np.abs(particle.pid)
    return (pid > 0 and pid < 7)

def IsLepton(particle):
    pid = np.abs(particle.pid)
    return (pid > 10 and pid < 19)

def IsGluon(particle):
    pid = np.abs(particle.pid)
    return (pid in [9,21])

def SelectFinalState(event):
    selection = (npyth.STATUS == 1) & (npyth.ABS_PDG_ID > -1)
    return event.all(selection=selection, return_hepmc=False)

def GatherQuarks(particle, particle_list=None):
    if(particle_list is None): particle_list = []
    if(IsQuark(particle)):
        particle_list.append(particle)
    elif(not IsGluon(particle)):
        children = particle.children(return_hepmc=True)
        if(len(children) != 0):
            for child in children:
                particle_list = GatherQuarks(child,particle_list)
    return particle_list

def SelectSimplestHadronic(event,truth_selection):
    truth_arr = [event.first(selection=x, return_hepmc=True) for x in truth_selection]
    particles = []
    for par in truth_arr:
        particles += GatherQuarks(par)
    return particles # TODO: Convert from HepMC to numpy
