import numpy as np
import pyhepmc as pyhep
from util.hepmc import HepMCArrayToNumpy

# ==============================
# These are particle selectors, that the user should
# use (via the configuration). Some are built to
# use specialized algorithms, which are contained
# in selection_algos.py.
# ==============================

# Called on a PythiaWrapper (that has generated an event),
# this class returns the index of the first particle in the event listing
# that has a matching PDG ID code and status.
class FirstSelector:
    def __init__(self,status,pdgid,hadronization=True):
        self.SetStatus(status)
        self.SetPdgId(pdgid)
        self.hadronization = hadronization
        self.selection_status = True

    def SetStatus(self,status):
        self.status = status

    def SetPdgId(self,pdgid):
        self.pdgid = pdgid

    def GetN(self):
        return 1

    def GetSelectionStatus(self):
        return self.selection_status

    def IsFixedLength(self):
        return True

    def SetHadronization(self,hadronization=True):
        self.hadronization = hadronization
        if(not hadronization):
            if np.abs(self.pdgid) in [1,2,3,4,5]:
                self.status = 1

    def __call__(self,pythia_wrapper):
        pdgid = pythia_wrapper.GetPdgId()
        status = pythia_wrapper.GetStatus(hepmc=True)
        self.selection_status = True

        pdgid_idx  = np.where(pdgid == self.pdgid)[0]
        status_idx = np.where(status == self.status)[0]

        intersection = np.intersect1d(pdgid_idx, status_idx)
        del pdgid_idx, status_idx
        if(len(intersection) == 0):
            self.selection_status = False
            return None
        return intersection[0]

    def Print(self):
        print('FirstSelector: status = {}, pdgid = {}'.format(self.status, self.pdgid))
        return

# This is a basic particle selection class -- it takes a list of FirstSelectors,
# each of which will select one particle. Thus it returns a list of selected particle indices,
# which should be of a fixed length.
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

    def __call__(self,pythia_wrapper):
        self.selection_status = True
        particle_list = []
        for x in self.selection_list:
            x.SetHadronization(self.hadronization)
            particle = x(pythia_wrapper)
            if(particle is not None): particle_list.append(particle)
            else: self.selection_status = False
        return np.sort(np.array(particle_list,dtype=int))

# This is a more advanced particle selection class -- it takes a particle selection algorithm,
# and returns that algorithm's output (truncated to some requested length n).
class AlgoSelection():
    def __init__(self,selection_algo,n, fixed_length=False):
        self.particle_selection_algo = selection_algo
        self.n = n
        self.fixed_length = fixed_length # whether or not this selector will always return the exact same number of particles
        self.selection_status = True

    def GetN(self):
        return self.n

    def GetSelectionStatus(self):
        return self.selection_status

    def IsFixedLength(self):
        return self.fixed_length

    def __call__(self,pythia_wrapper):
        self.selection_status, particle_list = self.particle_selection_algo(pythia_wrapper)
        if(len(particle_list) > self.n): particle_list = particle_list[:self.n]
        return np.sort(particle_list)

# This is a simple class for passing a list comprised of the above selectors.
# This allows one to construct more complex particle selections, implementing multiple algorithms.
class MultiSelection:
    def __init__(self,particle_selection_list):
        self.particle_selection_list = particle_selection_list
        self.n = np.sum([x.GetN() for x in particle_selection_list],dtype=int)

        self.fixed_length = True
        fixed_lengths = np.array([x.IsFixedLength() for x in self.particle_selection_list])
        if(False in fixed_lengths): self.fixed_length = False

        self.selection_status = True

    def GetN(self):
        return self.n

    def IsFixedLength(self):
        return self.fixed_length

    def GetSelectionStatus(self):
        return self.selection_status

    def __call__(self,pythia_wrapper):
        self.selection_status = True
        particle_lists = []
        statuses = []
        for i,selector in enumerate(self.particle_selection_list):
            particle_lists.append(selector(pythia_wrapper))
            statuses.append(selector.GetSelectionStatus())

        if(False in statuses): self.selection_status = False
        particle_list = np.unique(np.hstack(particle_lists))
        return np.sort(particle_list)

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