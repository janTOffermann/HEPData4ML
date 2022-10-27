import numpy as np
import pyhepmc as pyhep
from util.hepmc import HepMCArrayToNumpy

# Called on a PythiaWrapper (that has generated an event),
# this class returns the first particle in the event listing
# that has a matching PDG ID code and status, together with a list
# of daughter indices.
class FirstSelector:
    def __init__(self,status,pdgid,hadronization=True):
        self.SetStatus(status)
        self.SetPdgId(pdgid)

    def SetStatus(self,status):
        self.status = status

    def SetPdgId(self,pdgid):
        self.pdgid = pdgid

    def SetHadronization(self,hadronization=True):
        self.hadronization = hadronization
        if(not hadronization):
            if np.abs(self.pdgid) in [1,2,3,4,5]:
                self.status = 1

    def __call__(self,pythia_wrapper,return_hepmc=False):

        pdgid = pythia_wrapper.GetPdgId()
        status = pythia_wrapper.GetStatus(hepmc=True)

        pdgid_idx  = np.where(pdgid == self.pdgid)[0]
        status_idx = np.where(status == self.status)[0]

        intersection = np.intersect1d(pdgid_idx, status_idx)
        del pdgid_idx, status_idx
        if(len(intersection) == 0):
            return None

        idx = intersection[0]
        del intersection

        pmu = pythia_wrapper.GetPxPyPzE(indices=idx)
        pdgid = pdgid[idx]
        status = status[idx]
        daughters = pythia_wrapper.GetDaughtersSingle(index=idx)

        if(return_hepmc): # return a pyhepmc GenParticle
            v = pyhep.FourVector(*pmu) # px, py, pz, E
            particle = pyhep.GenParticle(v,pdgid,status)
            return particle

        else: # return a structured numpy array (E, px, py, pz, pdgid, status)
            particle = np.array(
                [(*pmu[1:],pmu[0],pdgid,status)],
                dtype = [('px','<f8'),('py','<f8'),('pz','<f8'),('e','<f8'),('pid','<i4'),('status','<i4')]
            )
            return particle

# This is a basic particle selection class -- it takes a list of FirstSelectors, each of which will select one particle.
class BasicSelection:
    def __init__(self,selection_list,hadronization=True):
        self.SetHadronization(hadronization)
        self.selection_list = selection_list
        self.n = len(self.selection_list)
        self.fixed_length = True # By design, this selector should always return a list of truth-level particles of some fixed length.
        self.selection_status = True

    def SetHadronization(self,hadronization=True):
        self.hadronization=hadronization

    def GetN(self):
        return self.n

    def GetSelectionStatus(self):
        return self.selection_status

    def IsFixedLength(self):
        return self.fixed_length

    def __call__(self,pythia_wrapper,return_hepmc=False):
        self.selection_status = True
        particle_list = []
        for x in self.selection_list:
            x.SetHadronization(self.hadronization)
            particle = x(pythia_wrapper,return_hepmc=return_hepmc)
            if(particle is not None): particle_list.append(particle)
            else: self.selection_status = False
        if(not return_hepmc): particle_list = np.concatenate(particle_list)
        return particle_list

# This is a more advanced particle selection class -- it takes a particle selection algorithm,
# and returns that algorithm's output (truncated to some requested length n).
class AlgoSelection():
    def __init__(self,selection_algo,n):
        self.particle_selection_algo = selection_algo
        self.n = n
        self.fixed_length = False

    def GetN(self):
        return self.n

    def IsFixedLength(self):
        return self.fixed_length

    def __call__(self,pythia_wrapper,return_hepmc=False):
        status, particle_list = self.particle_selection_algo(pythia_wrapper)
        if(len(particle_list) > self.n): particle_list = particle_list[:self.n]
        return particle_list


# selections = {
#     't->Wb': [
#         BasicSelection(22,6), # top from ttbar
#         BasicSelection(23,5), # bottom from t -> W b
#         BasicSelection(22,24) # W boson from t -> W b
#     ],
#     't->Wb_nohad': [ # TODO: Should be redundant thanks to hadronization toggle
#         BasicSelection(22,6), # top from ttbar
#         BasicSelection(1,5), # bottom from t -> W b
#         BasicSelection(22,24) # W boson from t -> W b
#     ],
#     'Wb': [
#         BasicSelection(23,5), # bottom from t -> W b
#         BasicSelection(22,24) # W boson from t -> W b
#     ],
#     'Wb_nohad': [ # TODO: Should be redundant thanks to hadronization toggle
#         BasicSelection(1,5), # bottom from t -> W b
#         BasicSelection(22,24) # W boson from t -> W b
#     ],
#     'W': [
#         BasicSelection(22,24) # W boson from t -> W b
#     ],
#     'b': [
#         BasicSelection(23,5), # bottom from bbar production
#     ],
#     'bbar': [
#         BasicSelection(23,-5), # anti-bottom from bbar production
#     ]
# }

# # This is a basic truth particle selection class. It takes an entry from the above dictionary.
# class TruthSelection:
#     def __init__(self,selection,hadronization=True):
#         self.SetHadronization(hadronization)
#         self.selection = selection
#         self.n = len(selection)
#         self.fixed_length = True # By design, this selector should always return a list of truth-level particles of some fixed length.

#     def SetHadronization(self,hadronization=True):
#         self.hadronization=hadronization

#     def GetN(self):
#         return self.n

#     def IsFixedLength(self):
#         return self.fixed_length

#     def __call__(self,event,return_hepmc=False):
#         selection_list = selections[self.selection]
#         particle_list = []
#         for x in selection_list:
#             x.SetHadronization(self.hadronization)
#             particle_list.append(x(event,return_hepmc=return_hepmc))
#         if(not return_hepmc): particle_list = np.concatenate(particle_list, axis=0)
#         # if(not return_hepmc): particle_list = HepMCArrayToNumpy(particle_list)
#         return particle_list

# This is a wrapper for particle selection methods from particle_selection.
# It truncates the output list to some fixed length n, so there will be at most n particles.
# (This does not perform any kind of zero-padding).
class AdvTruthSelection:
    def __init__(self,particle_selection, n):
        self.particle_selection = particle_selection
        self.n = n
        self.fixed_length = False # By design, this selector returns a variable-length list of truth-level particles (with some maximum length).

    def GetN(self):
        return self.n

    def IsFixedLength(self):
        return self.fixed_length

    def __call__(self,event,return_hepmc=False):
        self.particle_selection.SetHepMC(return_hepmc)
        status, particle_list = self.particle_selection(event)

        # Optional truncation.
        if(len(particle_list) > self.n):
            particle_list = particle_list[:self.n]

        # print('\nParticle list: {}'.format(len(particle_list)))
        # print(particle_list)
        return particle_list

# This class should be passed a list of selectors (TruthSelections, and/or AdvTruthSelections).
# It acts as the union (combination) of all the selectors in this list.
class MultiTruthSelection:
    def __init__(self,particle_selection_list):
        self.particle_selection_list = particle_selection_list
        self.n = np.sum([x.GetN() for x in particle_selection_list],dtype=int)
        self.fixed_length = False

    def GetN(self):
        return self.n

    def IsFixedLength(self):
        return self.fixed_length

    def __call__(self,event,return_hepmc=False):
        particle_lists = []
        for i,selector in enumerate(self.particle_selection_list):
            particle_lists.append(selector(event,return_hepmc))

        particle_list = np.hstack(particle_lists)
        return particle_list