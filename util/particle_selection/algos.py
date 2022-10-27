import numpythia as npyth # for defining selections
import numpy as np
from util.hepmc import HepMCArrayToNumpy

def IsQuark(particle, use_hepmc=True):
    if(use_hepmc): pid = np.abs(particle.pid)
    else: pid = np.abs(particle[-2]) # numpy
    return (pid > 0 and pid < 7)

def IsLepton(particle, use_hepmc=True):
    if(use_hepmc): pid = np.abs(particle.pid)
    else: pid = np.abs(particle[-2]) # numpy
    return (pid > 10 and pid < 19)

def IsGluon(particle, use_hepmc=True):
    if(use_hepmc): pid = np.abs(particle.pid)
    else: pid = np.abs(particle[-2]) # numpy
    return (pid in [9,21])

def IsPhoton(particle, use_hepmc=True):
    if(use_hepmc): pid = np.abs(particle.pid)
    else: pid = np.abs(particle[-2]) # numpy
    return (pid == 22)

def IsNeutrino(particle, use_hepmc=True):
    if(use_hepmc): pid = np.abs(particle.pid)
    else: pid = np.abs(particle[-2]) # numpy
    return (pid in [12, 14, 16, 18])

def IsBoson(particle, use_hepmc=True):
    if(use_hepmc): pid = np.abs(particle.pid)
    else: pid = np.abs(particle[-2]) # numpy
    return (pid in [23, 24, 25])

def IsStable(particle, use_hepmc=True):
    if(use_hepmc): status = np.abs(particle.status)
    else: status = particle[-1] # numpy
    return (status == 1)

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

def GatherStableDaughters(particle):
    particle_list = []
    if(particle.status == 1):
        particle_list.append(particle)
    else:
        children = particle.children(return_hepmc=True)
        if(len(children) != 0):
            for child in children:
                particle_list += GatherStableDaughters(child)
    return particle_list