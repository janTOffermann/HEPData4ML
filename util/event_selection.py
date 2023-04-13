import numpy as np
from util.calcs import DeltaR2Vectorized

# Keep final-state particles within a certain distance of the selected truth particles.
# We should be generous with our selection to avoid throwing out things that
# have any chance of making their way into our jets later on. This is determined
# by the distance parameter.
class TruthDistanceSelection:
    def __init__(self,distance=2.4, n_truth=None):
        self.SetDistance(distance)
        self.SetNTruth(n_truth)

    def SetDistance(self,distance):
        self.distance = distance
        self.d2 = np.square(self.distance)

    def SetNTruth(self,n):
        self.n = n

    def __call__(self,pythia_wrapper,final_state_indices,truth_indices):
        if(self.n is None): self.n = len(truth_indices)
        fs_particles    = np.array([[pythia_wrapper.GetTheta(x), pythia_wrapper.GetPhi(x)] for x in final_state_indices   ],dtype=np.dtype('f8'))
        truth_particles = np.array([[pythia_wrapper.GetTheta(x), pythia_wrapper.GetPhi(x)] for x in truth_indices[:self.n]],dtype=np.dtype('f8'))

        distances = DeltaR2Vectorized(fs_particles, truth_particles)
        min_distances = np.min(distances,axis=1)
        selected = (np.abs(min_distances) < np.abs(self.d2))
        return final_state_indices[selected]
