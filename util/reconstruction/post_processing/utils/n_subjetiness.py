import itertools
import ROOT as rt
import h5py as h5
import numpy as np
from typing import TYPE_CHECKING, Optional, Union
from util.reconstruction.post_processing.utils.postprocessor_base import PostProcessorBase

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

# TODO: Something seems broken, ghost is (occasionally?) throwing off <jet_name>.Constituents.Collection.
#       Appears at end of Constituents.Pmu, but not end of that branch.

class NSubjetiness(PostProcessorBase):
    """
    Compute the N-subjetiness of a jet.

    See: https://arxiv.org/abs/1011.2268 [JHEP 03 (2011) 015]
    """

    def __init__(self,N=2):
        super().__init__()

        self.N = N
        self.jetdef = None
        self.name = 'NSubjetiness'
        self.branch_name = None
        self.print_prefix = '{}'.format(self.name)

        self.subjet_vecs = None
        self.rvec = rt.Math.PtEtaPhiMVector()
        self._allocate_buffer()

        self.citations = {
            "NSubjetiness":
            """
@article{Thaler:2010tr,
    author = "Thaler, Jesse and Van Tilburg, Ken",
    title = "{Identifying Boosted Objects with N-subjettiness}",
    eprint = "1011.2268",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "MIT-CTP-4191",
    doi = "10.1007/JHEP03(2011)015",
    journal = "JHEP",
    volume = "03",
    pages = "015",
    year = "2011"
}
            """
        }

    def ModifyInitialization(self,obj:'JetFinder'):
        """
        This function will modify the initialization so that
        the ghost vectors are generated and loaded into memory.
        """
        import fastjet as fj # can do this since JetFinder ensures this is set up and cached correctly

        # NOTE: obj_name_input effectively unused, but want to pass info through
        if(self.obj_name_input is None):
            self.SetInputObjectName(obj.jet_name)
        self.obj_name_output = self.obj_name_input # doesn't modify jet -> pass through object name

        self.jetdef = fj.JetDefinition(fj.kt_algorithm, obj.radius) # TODO: Same radius?

    def _allocate_buffer(self):
        self.subjet_vecs = [
            rt.Math.PxPyPzEVector(0.,0.,0.,0.)
            for j in range(self.N)
        ]
        return

    def _compute_tau(self,obj : 'JetFinder',i,subjets):
        t = 0.

        # fetch the pt of the jet constituents
        pt = obj.constituent_vectors_cyl[i][:,0]
        d0 = np.sum(pt) * obj.radius

        for j,jet in enumerate(subjets):
           self.subjet_vecs[j].SetCoordinates(jet.px(),jet.py(),jet.pz(),jet.E())

        for j,constituent_vec in enumerate(obj.constituent_vectors_cyl[i]):
            self.rvec.SetCoordinates(*constituent_vec) # pt, eta, phi, m
            dr = np.min([rt.Math.VectorUtil.DeltaR(self.rvec,subjet) for subjet in self.subjet_vecs]) # TODO: Is DeltaR2 less intensive than DeltaR? Might avoid some calls to square root.
            t += pt[j] * dr # NOTE: DeltaR here uses eta, not rapidity
        return t / d0

    def ModifyJets(self, obj : 'JetFinder'):
        """
        Compute tau_N for the jets. This involves re-clustering with the exclusive kt algorithm.
        """
        import fastjet as fj # NOTE: In practice, fastjet will have been initialized already by JetFinder. Can similarly do this in Softdrop

        # Fetch the jet constituents. This fills obj.constituent_indices -- good to do for safety.
        obj._fetchJetConstituents()

        # store tags in a dictionary, where the keys are the jet indices from obj
        self.tau = {key:0. for key in obj.jets_dict.keys()}

        for i,jet in obj.jets_dict.items():

            # Recluster the jet to exactly N subjets, using the exclusive kt algorithm.
            kt_cs = fj.ClusterSequence(jet.constituents(), self.jetdef) # member of class, otherwise goes out-of-scope when ref'd later
            subjets = kt_cs.exclusive_jets_up_to(self.N)
            # print('For [{}], type(subjets) = {}, subjets = {}'.format(i,type(subjets),subjets))

            # Compute tau_N
            self.tau[i] = self._compute_tau(obj,i,subjets)
        return

    def ModifyWrite(self,obj : 'JetFinder'):
        self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
        self._addValueToBuffer(obj)

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        """
        #NOTE: I considered saving the specific indices of the particles to which
        #      each jet is ghost-associated, but I think that gets a bit complicated
        #      and it's not clear that it would be worthwhile.
        if(self.branch_name is None):
            self.branch_name = '{}.{}.Tau{}'.format(obj.jet_name,self.name,self.N) # NOTE: use of obj.jet_name
        if(self.branch_name not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.branch_name,ndim=1,dtype=float)
        return

    def _addValueToBuffer(self,obj : 'JetFinder'):
        """
        Adds the tau_N value to the buffer, for writing.
        """
        obj.output_buffer.set(self.branch_name,obj._i,self.tau)
