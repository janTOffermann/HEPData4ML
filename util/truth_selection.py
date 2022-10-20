import numpythia as npyth # for defining selections
import numpy as np
import pyhepmc_ng as pyhep

class BasicSelection:
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

    def __call__(self,event,return_hepmc=False):
        selection = (npyth.STATUS == self.status) & (npyth.PDG_ID == self.pdgid)
        return event.first(selection=selection,return_hepmc=return_hepmc)

selections = {
    't->Wb': [
        BasicSelection(22,6), # top from ttbar
        BasicSelection(23,5), # bottom from t -> W b
        BasicSelection(22,24) # W boson from t -> W b
    ],
    't->Wb_nohad': [ # TODO: Should be redundant thanks to hadronization toggle
        BasicSelection(22,6), # top from ttbar
        BasicSelection(1,5), # bottom from t -> W b
        BasicSelection(22,24) # W boson from t -> W b
    ],
    'Wb': [
        BasicSelection(23,5), # bottom from t -> W b
        BasicSelection(22,24) # W boson from t -> W b
    ],
    'Wb_nohad': [ # TODO: Should be redundant thanks to hadronization toggle
        BasicSelection(1,5), # bottom from t -> W b
        BasicSelection(22,24) # W boson from t -> W b
    ],
    'b': [
        BasicSelection(23,5), # bottom from bbar production
    ],
    'bbar': [
        BasicSelection(23,-5), # anti-bottom from bbar production
    ],
    'W': [
        BasicSelection(22,24) # W boson from t -> W b
    ]
}

# This is a basic truth particle selection class. It takes an entry from the above dictionary.
class TruthSelection:
    def __init__(self,selection,hadronization=True):
        self.SetHadronization(hadronization)
        self.selection = selection

    def SetHadronization(self,hadronization=True):
        self.hadronization=hadronization

    def __call__(self,event,return_hepmc=False):
        selection_list = selections[self.selection]
        particle_list = []
        for x in selection_list:
            x.SetHadronization(self.hadronization)
            particle_list.append(x(event,return_hepmc=return_hepmc))
        if(not return_hepmc): particle_list = np.concatenate(particle_list, axis=0)
        return particle_list

# This is a wrapper for particle selection methods from particle_selection.
# It truncates the output list to some fixed length n.
class AdvTruthSelection:
    def __init__(self,particle_selection, n):
        self.particle_selection = particle_selection
        self.n = n

    def __call__(self,event,return_hepmc=False):
        self.particle_selection.SetHepMC(return_hepmc)
        status, particle_list = self.particle_selection(event)
        # Optional truncation.
        if(len(particle_list) > self.n):
            particle_list = particle_list[:self.n]

        print('\nParticle list: {}'.format(len(particle_list)))
        print(particle_list)
        return particle_list
