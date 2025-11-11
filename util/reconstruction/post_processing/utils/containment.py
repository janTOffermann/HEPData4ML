import numpy as np
import ROOT as rt
from typing import TYPE_CHECKING
from util.reconstruction.post_processing.utils.postprocessor_base import PostProcessorBase

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

class ContainmentTagger(PostProcessorBase):
    """
    Tags jets based on a deltaR check, against (user-specified) elements of some input collection specified by "key".
    For more advanced containment checks, consider GhostAssociation.
    """
    def __init__(self,key,indices,delta_r=None, mode='filter',use_rapidity=True, tag_name=None):
        super().__init__()

        self.key = key
        self.vec_key = '{}.Pmu_cyl'.format(key)
        self.indices = indices
        self.vecs = None # will store necessary vectors in memory #TODO: Make this more memory-friendly and read from file as needed?

        self.mode = mode
        assert self.mode in ['tag','filter']
        self.tag_name = tag_name

        self.SetRadius(delta_r)

        self.use_rapidity = use_rapidity

        self.tags = None

        self.name = 'ContainmentTagger'
        self.print_prefix = '\n\t{}'.format(self.name)

        # Transient, per-jet variables
        self.tag_status = False
        self.vec1 = rt.Math.PtEtaPhiMVector()
        self.vec2 = rt.Math.PtEtaPhiMVector()

    def SetRadius(self,dr:float):
        self.radius2 = np.square(dr)

    def _compute_distance2(self,v1,v2):
        self.vec1.SetCoordinates(*v1)
        self.vec2.SetCoordinates(*v2)

        if(self.use_rapidity):
            dphi = rt.Math.VectorUtil.DeltaPhi(self.vec1,self.vec2)
            dy = self.vec2.Rapidity() - self.vec1.Rapidity()
            return np.square(dphi) + np.square(dy)
        else:
            return rt.Math.VectorUtil.DeltaR2(self.vec1,self.vec2)

    def _tag(self,obj : 'JetFinder', pmu_cyl):
        status = True
        # Fetch the vectors from the input buffer; we don't assume that they
        # are necessarily in obj.input_collection_arrays (they likely are not).
        vecs = np.array(obj.input_buffer[self.vec_key])[self.indices] # reminder: using cylindrical

        if(vecs.ndim == 1):
            distance2 = self._compute_distance2(pmu_cyl,vecs)
            if(distance2 > self.radius2):
                status = False
        else:

            for i,vec in enumerate(vecs):
                distance2 = self._compute_distance2(pmu_cyl,vec)
                if(distance2 > self.radius2):
                    status = False
                    break

        self.tag_status = status
        return

    def ModifyInitialization(self,obj : 'JetFinder'):
        """
        This function will modify the initialization so that
        the necessary input vectors are loaded into memory, if
        not already present.
        """
        if(self.obj_name_input is None):
            self.SetInputObjectName(obj.jet_name)
        self.obj_name_output = self.obj_name_input # doesn't modify jet -> pass through object name

        if(self.indices is not None):
            self.indices = np.atleast_1d(self.indices)

        # Fetch the 4-vector key, and make sure its data is loaded.
        # Note the use of cylindrical coordinates!
        if(self.vec_key not in obj.input_collection_arrays.keys()):
            obj.input_buffer.read_branch(self.vec_key)

        if(self.mode=='tag'):
            self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will tag jets, and fill the corresponding branches.
        """

        pmu_cyl_key = '{}.Pmu_cyl'.format(self.obj_name_input)
        pmu_cyl_dict = obj.output_buffer.get(pmu_cyl_key,filter=obj.jet_ordering)
        self.tags = {i:False for i in pmu_cyl_dict.keys()}

        for key,pmu_cyl in obj.pmu_cyl_dict.items():
            self._tag(obj,pmu_cyl) # fills self.tag_status
            self.tags[key] = self.tag_status

        # For filtering mode, we modify obj.jet_ordering
        if(self.mode=='filter'):
            obj.jet_ordering = [key for key in obj.jet_ordering if self.tags[key]]
            # TODO: With the shift to how data I/O is handled, do we actually need to call the below funcs anymore?
            # obj._updateJetDictionary()
            # obj._jetsToVectors()
            # obj._fetchJetConstituents()

    def ModifyWrite(self,obj : 'JetFinder'):
        if(self.mode=='filter'):
            return # do nothing
        else:
            self._addFlagToBuffer(obj)

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Used if self.mode=='tag', in which case we're writing all jets,
        and including a new branch to indicate whether or not a jet is
        JH-tagged.
        """
        self._createBranchNames(obj)

        if(self.tag_name not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.tag_name,ndim=1,dtype=bool)
        return

    def _createBranchNames(self,obj : 'JetFinder'):
        if(self.tag_name is None):
            self.tag_name = '{}.{}'.format(self.obj_name_output,self.name)

    def _addFlagToBuffer(self,obj : 'JetFinder'):
        """
        Adds the containment tags to the buffer, for writing.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.output_buffer.set(self.tag_name,obj._i,self.tags)