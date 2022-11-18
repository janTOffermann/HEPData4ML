import numpy as np
from util.calcs import DeltaR2

class GetNearestJet:
    def __init__(self,truth_code=6, max_dr=0.8, use_hepmc=True):
        self.SetTruthCode(truth_code)
        self.SetMaxDeltaR(max_dr)
        self.SetUseHepMC(use_hepmc)

    def SetTruthCode(self,code):
        self.truth_code = code

    def SetMaxDeltaR(self,dr):
        self.max_dr = dr

    def SetUseHepMC(self,flag):
        self.use_hepmc = flag

    def __call__(self,**kwargs):
        truth_particles = kwargs['truth']
        jets = kwargs['jets']

        # Get the truth-level particle for which we want to match a jet.
        if(self.use_hepmc): pdg_codes = np.array([x.pid for x in truth_particles])
        else: pdg_codes = truth_particles[:]['pdgid']
        if(self.truth_code not in pdg_codes):
            print('Error: No truth-level particle found with code {}.'.format(self.truth_code))
            return -1 # assert False

        truth_idx = np.where(pdg_codes == self.truth_code)[0][0]
        truth = truth_particles[truth_idx]

        if(self.use_hepmc):
            truth_eta = truth.momentum.eta()
            truth_phi = truth.momentum.phi()
        else:
            truth_eta = truth['eta']
            truth_phi = truth['phi']

        # Now find distance between the truth particle and each of the jets, and pick the closest jet.
        dr2 = np.array([DeltaR2(truth_eta,truth_phi,x.eta(),x.phi()) for x in jets]).flatten()

        # Optional check on the distance of the nearest jet.
        if(self.max_dr > 0.):
            min_dr2 = np.min(dr2)
            if(min_dr2 > self.max_dr * self.max_dr): return -1 # no jet within requested radius
        return np.argmin(dr2)

class GetLeadingJet:
    def __init__(self):
        pass
    def __call__(self,**kwargs):
        jets = kwargs['jets']
        pt = np.array([x.pt() for x in jets])
        return np.argmax(pt)