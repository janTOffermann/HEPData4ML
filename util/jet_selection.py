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

class GetNearestJetWithContainment:
    """
    Like GetNearestJet but in addition, this algorithm
    will require that the provided, truth-level particles,
    within the given index range, fall within the (given)
    "containment radius" of this jet.
    For example, with the appropriate choice of truth selection
    (for the truth particles) one could use this to select jets
    nearest a top quark in (eta,phi), and then require that the
    stable decay products of the W from t->Wb are all contained
    within the jet radius.

    Note that this containment simply uses a distance check
    in (eta, phi). It does not explicitly check that each truth-
    level particle explicitly corresponds with a jet constituent.
    """
    def __init__(self,truth_code=6, max_dr=0.8, truth_containment_start_idx=0, truth_containment_stop_idx=-1,containment_radius=None,use_hepmc=True):
        self.jet_selector = GetNearestJet(truth_code,max_dr,use_hepmc) # this will first run the GetNearestJet algo, then apply an additional cut.
        self.truth_start_idx = truth_containment_start_idx
        self.truth_stop_idx = truth_containment_stop_idx
        self.containment_radius = containment_radius
        self.use_hepmc = use_hepmc
        if(containment_radius is None):
            self.containment_radius = max_dr

    def __call__(self,**kwargs):
        truth_particles = kwargs['truth']
        jets = kwargs['jets']
        jet_idx = self.jet_selector(**kwargs)
        if(jet_idx < 0): return -1 # GetNearestJet failed to find a jet
        jet = jets[jet_idx]

        # Now check the distance between the selected jet and each truth particle within our given indices
        n_truth = len(truth_particles)
        if(self.truth_stop_idx < 0): self.truth_stop_idx = len(truth_particles) + self.truth_stop_idx
        if(self.truth_stop_idx >= n_truth): self.truth_stop_idx = n_truth - 1
        truth_particles = truth_particles[self.truth_start_idx:self.truth_stop_idx+1]

        if(self.use_hepmc):
            truth_eta = [t.momentum.eta() for t in truth_particles]
            truth_phi = [t.momentum.phi() for t in truth_particles]
        else:
            truth_eta = [t['eta'] for t in truth_particles]
            truth_phi = [t['phi'] for t in truth_particles]

        dr2 = np.array([DeltaR2(jet.eta(),jet.phi(),truth_eta[i],truth_phi[i]) for i in range(len(truth_particles))]).flatten()

        if(np.max(dr2) > self.containment_radius * self.containment_radius): return -1
        return jet_idx

class GetLeadingJet:
    def __init__(self):
        pass
    def __call__(self,**kwargs):
        jets = kwargs['jets']
        pt = np.array([x.pt() for x in jets])
        return np.argmax(pt)