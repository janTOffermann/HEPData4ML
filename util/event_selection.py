import numpy as np
from util.calcs import DeltaR2Vectorized
from util.hepmc import RestructureParticleArray

# Keep final-state particles within a certain distance of the selected truth particles.
# We should be generous with our selection to avoid throwing out things that
# have any chance of making their way into our jets later on. This is determined
# by the alpha parameter, which should always be greater than 1.
class TruthDistanceSelection:
    def __init__(self,distance=2.4):
        self.SetDistance(distance)

    def SetDistance(self,distance):
        self.distance = distance

    def __call__(self,arr,arr_truth):
        arr_truth_tmp = RestructureParticleArray(arr_truth)
        arr_tmp = RestructureParticleArray(arr)
        # arr_truth = RestructureParticleArray(arr_truth)
        # arr = RestructureParticleArray(arr)
        dr2_limit = np.square(self.distance)

        distances = DeltaR2Vectorized(arr_tmp[:,-2:], arr_truth_tmp[:,-2:])
        min_distances = np.min(distances,axis=1)
        selected = (np.abs(min_distances) < np.abs(dr2_limit))
        arr = arr[selected]

        status = len(arr) > 0
        del arr_truth_tmp, arr_tmp, distances, min_distances, selected
        return status, arr,arr_truth