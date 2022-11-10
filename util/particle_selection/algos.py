# Some very generic/simple algorithms, which are
# used in selection_algos.py. Keeping them here
# to limit code clutter elsewhere.

import numpy as np

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

def IsQuark(pythia_wrapper=None, idx=None, p=None, hepmc_particle=None):
    if(hepmc_particle is not None):
        pid = np.abs(hepmc_particle.pid)
    else:
        if(p is None): p = pythia_wrapper.GetParticle(idx)
        pid = np.abs(p[-2])
    return (pid > 0 and pid < 7)

def IsLepton(pythia_wrapper=None, idx=None, p=None, hepmc_particle=None):
    if(hepmc_particle is not None):
        pid = np.abs(hepmc_particle.pid)
    else:
        if(p is None): p = pythia_wrapper.GetParticle(idx)
        pid = np.abs(p[-2])
    return (pid > 10 and pid < 19)

def IsGluon(pythia_wrapper=None, idx=None, p=None, hepmc_particle=None):
    if(hepmc_particle is not None):
        pid = np.abs(hepmc_particle.pid)
    else:
        if(p is None): p = pythia_wrapper.GetParticle(idx)
        pid = np.abs(p[-2])
    return (pid in [9,21])

def IsPhoton(pythia_wrapper=None, idx=None, p=None, hepmc_particle=None):
    if(hepmc_particle is not None):
        pid = np.abs(hepmc_particle.pid)
    else:
        if(p is None): p = pythia_wrapper.GetParticle(idx)
        pid = np.abs(p[-2])
    return (pid == 22)

def IsNeutrino(pythia_wrapper=None, idx=None, p=None, hepmc_particle=None):
    if(hepmc_particle is not None):
        pid = np.abs(hepmc_particle.pid)
    else:
        if(p is None): p = pythia_wrapper.GetParticle(idx)
        pid = np.abs(p[-2])
    return (pid in [12, 14, 16, 18])

def IsBoson(pythia_wrapper=None, idx=None, p=None, hepmc_particle=None):
    if(hepmc_particle is not None):
        pid = np.abs(hepmc_particle.pid)
    else:
        if(p is None): p = pythia_wrapper.GetParticle(idx)
        pid = np.abs(p[-2])
    return (pid in [23, 24, 25])

def IsStable(pythia_wrapper=None, idx=None, p=None, hepmc_particle=None):
    if(hepmc_particle is not None):
        status = np.abs(hepmc_particle.status)
    else:
        if(p is None): p = pythia_wrapper.GetParticle(idx)
        status = p[-1]
    return (status == 1)

# # Given the index of a particle in the event listing
# # of a Pythia wrapper, traverse down the event listing
# # and gather the nearest quark daughters.
class GatherQuarks:
    def __init__(self):
        self.particle_list = []
        # self.searched_particle_list = [] # to avoid loops -- though this was probably due to a bug during early development!

    def __run(self,pythia_wrapper,idx):
        # self.searched_particle_list.append(idx)
        daughters = pythia_wrapper.GetDaughtersSingle(idx)

        # Determine which daughters are quarks -- these get added to particle_list.
        # Also determine which ones are not -- we traverse down these parts of the tree.
        quark_indices = np.array([IsQuark(pythia_wrapper,x) for x in daughters],dtype=bool)
        not_quark_indices = np.array([not x for x in quark_indices],dtype=bool)

        quark_daughters = daughters[quark_indices]
        not_quark_daughters = daughters[not_quark_indices]

        for q in quark_daughters:
            if(q not in self.particle_list):
                self.particle_list.append(q)

        for nq in not_quark_daughters:
            # if(nq in self.searched_particle_list): continue
            self.__run(pythia_wrapper,nq)

    def __call__(self,pythia_wrapper,idx):
        self.__run(pythia_wrapper,idx)
        return np.array(self.particle_list,dtype=int)

class GatherStableDaughters:
    def __init__(self):
        self.particle_list = []
        # self.searched_particle_list = []

    def __run(self,pythia_wrapper,idx):
        # self.searched_particle_list.append(idx)
        daughters = pythia_wrapper.GetDaughtersSingle(idx)

        # Determine which daughters are stable.
        # Also determine which ones are not.
        stable_indices = np.array([IsStable(pythia_wrapper,x) for x in daughters],dtype=bool)
        not_stable_indices = np.array([not x for x in stable_indices],dtype=bool)

        stable_daughters = daughters[stable_indices]
        not_stable_daughters = daughters[not_stable_indices]

        for d in stable_daughters:
            if(d not in self.particle_list):
                self.particle_list.append(d)

        for nd in not_stable_daughters:
            # if(nd in self.searched_particle_list): continue
            self.__run(pythia_wrapper,nd)

    def __call__(self,pythia_wrapper,idx):
        self.__run(pythia_wrapper,idx)
        return np.array(self.particle_list,dtype=int)