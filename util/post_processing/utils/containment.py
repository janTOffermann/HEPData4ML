import os,glob,pathlib
import subprocess as sub
import numpy as np
import ROOT as rt
import h5py as h5
from util.calcs import embed_array_inplace
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.post_processing.jets import JetFinder


class ContainmentTagger:
    """
    Tags jets based on a deltaR check, against (user-specified) elements of some input collection specified by "key".
    For more advanced containment checks, consider GhostAssociation.
    """
    def __init__(self,key,indices,delta_r=None, mode='filter',use_rapidity=True, tag_name=None):

        self.key = key
        self.vec_key = '{}.Pmu_cyl'.format(key)
        self.indices = indices
        self.vecs = None # will store necessary vectors in memory #TODO: Make this more memory-friendly and read from file as needed?

        self.mode = mode
        assert self.mode in ['tag','filter']
        self.tag_name = tag_name

        # Transient, per-jet variables
        self.tag_status = False

        self.SetRadius(delta_r)

        self.use_rapidity = use_rapidity

        self.tags = None

        self.print_prefix = '\n\t\tContainmentTagger'
        self.citations = {}

    def GetCitations(self):
        return self.citations

    def SetRadius(self,dr:float):
        self.radius2 = np.square(dr)

    def _compute_distance2(self,v1,v2):
        vec1 = rt.Math.PtEtaPhiMVector(*v1)
        vec2 = rt.Math.PtEtaPhiMVector(*v2)

        if(self.use_rapidity):
            dphi = rt.Math.VectorUtil.DeltaPhi(vec1,vec2)
            dy = vec2.Rapidity() - vec1.Rapidity()
            return np.square(dphi) + np.square(dy)
        else:
            return rt.Math.VectorUtil.DeltaR2(vec1,vec2)

    def _tag(self,obj : 'JetFinder', key:int):
        status = True

        vecs = obj.input_collection_arrays[self.vec_key][obj._i] # reminder: using cylindrical
        jet_vec = obj.jet_vectors_cyl[key]

        # Compute distances. Deal with cases of "vecs" being multiple vectors, or a single one.

        if(vecs.ndim == 1):
            distance2 = self._compute_distance2(jet_vec,vecs)
            if(distance2 > self.radius2):
                status = False
        else:

            for i,vec in enumerate(vecs):
                distance2 = self._compute_distance2(jet_vec,vec)
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
        if(self.indices is not None):
            self.indices = np.atleast_1d(self.indices)

        # Fetch the 4-vector key, and make sure its data is loaded.
        # Note the use of cylindrical coordinates!
        if(self.vec_key not in obj.input_collection_arrays.keys()):
            # Read the necessary keys in.
            f = h5.File(obj.h5_file,'r')

            # Cylindrical coordinates four-momenta
            if(self.indices is not None):

                data = f[self.vec_key][:]
                if(data.ndim == 2):
                    self._print('Warning: indices != None, but branch {} has ndim = {}. Setting indices -> None.'.format(self.vec_key,data.ndim))
                    self.indices = None
                else:
                    vecs = f[self.vec_key][:][:,self.indices] # only loads the data needed -- indexing already done here
            if(self.indices is None):
                vecs = f[self.vec_key][:]

            obj.input_collection_arrays[self.vec_key] = vecs # NOTE: This won't affect jet clustering, as we don't modify obj.input_collections
            f.close()

    def ModifyInputs(self,obj : 'JetFinder'):
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will tag jets, and fill the corresponding branches.
        """
        import fastjet as fj # NOTE: In practice, fastjet will have been initialized already by JetFinder. Can similarly do this in Softdrop

        # Fetch the jet constituents, just to be safe -- this makes sure that they are up-to-date.
        obj._fetchJetConstituents()

        self.tags = {i:False for i in obj.jets_dict.keys()}

        for key in obj.jets_dict.keys():
            self._tag(obj,key) # fills self.tag_status
            self.tags[key] = self.tag_status

        if(self.mode=='filter'):
            obj.jet_ordering = [key for key in obj.jet_ordering if self.tags[key]]
            obj._updateJetDictionary()
            # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
            obj._jetsToVectors()
            obj._fetchJetConstituents()

    def ModifyConstituents(self, obj : 'JetFinder'):
        return

    def ModifyWrite(self,obj : 'JetFinder'):
        if(self.mode=='filter'):
            return # do nothing
        else:
            self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
            self._addFlagToBuffer(obj)

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Used if self.mode=='tag', in which case we're writing all jets,
        and including a new branch to indicate whether or not a jet is
        JH-tagged.
        """
        self._createBranchNames(obj)

        if(self.tag_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.tag_name,(obj.n_jets_max,),dtype=bool)
        return

    def _createBranchNames(self,obj : 'JetFinder'):
        if(self.tag_name is None):
            self.tag_name = '{}.ContainmentTagger'.format(obj.jet_name)

    def _addFlagToBuffer(self,obj : 'JetFinder'):
        """
        Adds the containment tags to the buffer, for writing.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.buffer.set(self.tag_name,obj._i,[self.tags[i] for i in obj.jet_ordering])

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return
