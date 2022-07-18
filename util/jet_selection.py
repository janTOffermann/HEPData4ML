import numpy as np
from util.calcs import DeltaR2

# Select the jet nearest to a truth-level top quark.
def GetTopJet(**kwargs):
    truth_particles = kwargs['truth']
    jets = kwargs['jets']
    top_code = 6

    # Determine the format of truth_particles.
    use_hepmc = False
    if('use_hepmc' in kwargs.keys()): use_hepmc = kwargs['use_hepmc']

    # Maximum distance allowed.
    max_dr = 0.8
    if('max_dr' in kwargs.keys()): max_dr = kwargs['max_dr']

    # Get the truth top.
    if(use_hepmc): pdg_codes = np.array([x.pid for x in truth_particles])
    else: pdg_codes = truth_particles[:]['pdgid']
    if(top_code not in pdg_codes):
        print('Error: No truth top found.')
        assert(False)

    top_idx = np.where(pdg_codes == top_code)[0][0]
    top = truth_particles[top_idx]

    # Now find distance between the top and each of the jets, and pick the closest jet.
    if(use_hepmc): dr2 = np.array([DeltaR2(top.momentum.eta(),top.momentum.phi(),x.eta(),x.phi()) for x in jets]).flatten()
    else: dr2 = np.array([DeltaR2(top['eta'],top['phi'],x.eta(),x.phi()) for x in jets]).flatten()

    # Optional check on the distance of the nearest jet.
    if(max_dr > 0.):
        min_dr2 = np.min(dr2)
        if(min_dr2 > max_dr * max_dr): return -1 # no jet within requested radius
    return np.argmax(-dr2)

# Select the leading (highest pT) jet.
def GetLeadingJet(**kwargs):
    jets = kwargs['jets']
    pt = np.array([x.pt() for x in jets])
    return np.argmax(pt)
