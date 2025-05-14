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

def IsQuark(hepev,idx):
    pid = np.abs(hepev.particles[idx].pid)
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

def IsStable(hepev,idx):
    return hepev.particles[idx].status == 1

def GetDaughtersSingle(hepev,idx):
    return [x.id - 1 for x in hepev.particles[idx].children] # id gives 1-indexing # TODO: double-check this

def GetDaughtersRecursive(hepev,idx):
    daughters = GetDaughtersSingle(hepev,idx)
    # get the daughters' daughters
    for daughter in daughters:
        daughters.append(GetDaughtersRecursive(hepev,daughter))
    return daughters


# # Given the index of a particle in the event listing
# # of a Pythia wrapper, traverse down the event listing
# # and gather the nearest quark daughters.
class GatherQuarks:
    def __init__(self):
        self.particle_list = []

    def __run(self,hepev,idx):
        daughters = GetDaughtersSingle(hepev,idx) # indices of the daughters of particle @ idx

        # Determine which daughters are quarks -- these get added to particle_list.
        # Also determine which ones are not -- we traverse down these parts of the tree.
        quark_indices = np.array([IsQuark(hepev,x) for x in daughters],dtype=bool)
        not_quark_indices = np.array([not x for x in quark_indices],dtype=bool)

        quark_daughters = daughters[quark_indices]
        not_quark_daughters = daughters[not_quark_indices]

        for q in quark_daughters:
            if(q not in self.particle_list):
                self.particle_list.append(q)

        for nq in not_quark_daughters:
            self.__run(hepev,nq)

    def __call__(self,hepev,idx):
        self.__run(hepev,idx)
        return np.array(self.particle_list,dtype=int)

class GatherStableDaughters:
    def __init__(self):
        self.particle_list = []

    def __run(self,hepev,idx):
        daughters = GetDaughtersSingle(hepev,idx) # indices of the daughters of particle @ idx

        # Determine which daughters are stable.
        # Also determine which ones are not.
        stable_indices = np.array([IsStable(hepev,x) for x in daughters],dtype=bool)
        not_stable_indices = np.array([not x for x in stable_indices],dtype=bool)

        stable_daughters = daughters[stable_indices]
        not_stable_daughters = daughters[not_stable_indices]

        for d in stable_daughters:
            if(d not in self.particle_list):
                self.particle_list.append(d)

        for nd in not_stable_daughters:
            self.__run(hepev,nd)

    def __call__(self,hepev,idx):
        self.particle_list.clear()
        self.__run(hepev,idx)
        return np.unique(np.array(self.particle_list,dtype=int))

class GatherDaughters:
    """
    Given some particle, look for its daughters (from decay).
    """
    def __init__(self,recursive=False):
        self.particle_list = []
        self.recursive=recursive

    def __run(self,hepev,idx):
        daughters = GetDaughtersRecursive(hepev,idx)

        # Check for "self" in daughter list (a particle with same PDG code).
        # (we don't want to be thrown off by Pythia's re-listings, or brem, etc.)
        daughter_particles = [hepev.particles[x] for x in daughters]
        daughter_pid = [x.pid for x in daughter_particles]
        this_pid = hepev.particles[idx].pid
        if(this_pid in daughter_pid):
            new_idx = daughters[np.where(daughter_pid==this_pid)[0]]
            self.__run(hepev,new_idx)
            return

        for d in daughters:
            if(d not in self.particle_list):
                self.particle_list.append(d)

    def __call__(self,hepev,idx):
        self.particle_list.clear()
        self.__run(hepev,idx)
        return np.array(self.particle_list,dtype=int)