import numpy as np

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
        if(self.status is not None):
            status_idx = np.where(status == self.status)[0]
            intersection = np.intersect1d(pdgid_idx, status_idx)
            del status_idx
        else: intersection = np.copy(pdgid_idx)
        del pdgid_idx
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
        self.SetN(n)
        self.fixed_length = fixed_length # whether or not this selector will always return the exact same number of particles
        self.selection_status = True

    def GetN(self):
        return self.n

    def SetN(self,n):
        self.n = n

    def GetSelectionStatus(self):
        return self.selection_status

    def IsFixedLength(self):
        return self.fixed_length

    def __call__(self,pythia_wrapper):
        self.selection_status, particle_list = self.particle_selection_algo(pythia_wrapper)
        if(len(particle_list) > self.n and self.n > 0): particle_list = particle_list[:self.n]
        return particle_list

# This is a simple class for passing a list comprised of the above selectors.
# This allows one to construct more complex particle selections, implementing multiple algorithms.
class MultiSelection:
    def __init__(self,particle_selection_list, enforce_unique=False):
        self.particle_selection_list = particle_selection_list
        self.n = np.sum([x.GetN() for x in particle_selection_list],dtype=int)

        self.fixed_length = True
        fixed_lengths = np.array([x.IsFixedLength() for x in self.particle_selection_list])
        if(False in fixed_lengths): self.fixed_length = False

        self.enforce_unique = enforce_unique
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
        # statuses = []
        for i,selector in enumerate(self.particle_selection_list):
            particle_list = selector(pythia_wrapper)
            if(particle_list is None): particle_list = -1 # a selector failed -- the status should be reported as False
            if(type(particle_list) not in [list,np.ndarray]): particle_list = np.array(particle_list,dtype=int)
            particle_lists.append(particle_list)
            individual_status = selector.GetSelectionStatus()
            if(not individual_status):
                self.selection_status = False
                break
            # statuses.append(selector.GetSelectionStatus())

        # if(False in statuses): self.selection_status = False
        # return np.hstack(particle_list).flatten()
        particle_list = np.hstack(particle_lists)
        if(self.enforce_unique): particle_list = np.unique(particle_list) # TODO: May have unexpected consequences for particle ordering
        return particle_list