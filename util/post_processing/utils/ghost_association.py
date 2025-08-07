import itertools
import ROOT as rt
import h5py as h5
import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.post_processing.jets import JetFinder

# TODO: Something seems broken, ghost is (occasionally?) throwing off <jet_name>.Constituents.Collection.
#       Appears at end of Constituents.Pmu, but not end of that branch.

class GhostAssociator():
    """
    Ghost-associates jets with particles labeled by `key`,
    at the given particle indices (within each event).

    See: https://arxiv.org/abs/0802.1188 [JHEP 04 (2008) 005]
    """

    def __init__(self,key,indices,mode='filter',tag_name=None):

        self.key = key
        self.vec_key = '{}.Pmu_cyl'.format(key)
        self.ghost_key = '{}.Ghost.Pmu'.format(key)
        self.indices = indices
        self.mode = mode
        assert self.mode in ['tag','filter']

        self.tag_name=tag_name
        self.tags = None
        self.print_prefix = '\n\t\tGhostAssociator'

    def _makeGhosts(self,vecs, a=1.0e-10):

        result = np.zeros(vecs.shape)
        for i,vec in enumerate(vecs):
            v = rt.Math.PtEtaPhiMVector(a, vec[1], vec[2], a)
            result[i] = np.array([v.E(),v.Px(),v.Py(),v.Pz()])
        return result

    def ModifyInitialization(self,obj : 'JetFinder'):
        """
        This function will modify the initialization so that
        the ghost vectors are generated and loaded into memory.
        """
        # self._initialize_fastjet()
        self.indices = np.atleast_1d(self.indices)

        # Fetch the 4-vector key, and make sure its data is loaded.
        # Note the use of cylindrical coordinates!
        if(self.vec_key not in obj.input_collection_arrays.keys()):
            # Read the necessary keys in.
            f = h5.File(obj.h5_file,'r')

            # Cylindrical coordinates four-momenta
            vecs = f[self.vec_key][:][:,self.indices] # only loads the data needed -- indexing already done here

            obj.input_collection_arrays[self.vec_key] = vecs
            f.close()

        # Now convert these to ghosts: send pT and m -> 0.
        # These are in Cartesian (E,px,py,pz)
        ghost_vecs = np.array([self._makeGhosts(v) for v in obj.input_collection_arrays[self.vec_key]])

        obj.input_collection_arrays[self.ghost_key] = ghost_vecs

        # if(self.mode=='tag'):
        #     self._initializeBuffer(obj)

    def ModifyInputs(self,obj : 'JetFinder'):
        """
        Puts the ghost vectors corresponding with `key`
        on the bottom of obj.input_vecs.
        """
        # Putting the ghosts on the bottom to avoid causing issues with fasjet.Pseudojet.user_index().
        ghost_vecs = obj.input_collection_arrays[self.ghost_key][obj._i] # ith event

        original_input_length = len(obj.input_vecs)

        obj.input_vecs = np.vstack([obj.input_vecs,ghost_vecs])

        for i in range(len(obj.input_vecs)):
            ghost_dict = {'GhostAssociation:Ghost':(i >= original_input_length)}
            obj.AddUserInfo(i,ghost_dict) # TODO: Fetch existing UserInfo first, and add this instead? Don't want to accidentally overwrite something else.
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will do one of two things, depending on self.mode:
        - self.mode == 'filter': Drop jets that aren't ghost-associated.
        - self.mode =='tag': Write outputs indicating whether or not the jet is ghost-associated.
        This function will also make sure that ghosts in the jets' constituent lists are removed.
        """
        import fastjet as fj # NOTE: In practice, fastjet will have been initialized already by JetFinder. Can similarly do this in Softdrop

        # Fetch the jet constituents. This fills obj.constituent_indices -- good to do for safety.
        obj._fetchJetConstituents()

        # store tags in a dictionary, where the keys are the jet indices from obj
        self.tags = {key:False for key in obj.jets_dict.keys()}

        for i,jet in obj.jets_dict.items():
            ghost_mask = [False]
            if(jet.has_constituents()):
                ghost_mask = np.array([x.python_info()['GhostAssociation:Ghost'] for x in jet.constituents()])
            self.tags[i] = np.sum(ghost_mask) > 0

            # For the jets with ghosts, modify them to remove the ghost --
            # we don't want to pass it to any further steps.
            if(self.tags[i]):
                obj.jets_dict[i] = fj.join([x for x in list(itertools.compress(list(jet.constituents()),~ghost_mask))])

        if(self.mode=='filter'):
            obj.jet_ordering = [key for key in obj.jet_ordering if self.tags[key]]
            obj._updateJetDictionary()
            # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
            obj._jetsToVectors()
            obj._fetchJetConstituents()

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
        ghost-associated.
        """
        #NOTE: I considered saving the specific indices of the particles to which
        #      each jet is ghost-associated, but I think that gets a bit complicated
        #      and it's not clear that it would be worthwhile.
        if(self.tag_name is None):
            self.tag_name = '{}.{}.GhostAssociated'.format(obj.jet_name,self.key)
        if(self.tag_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.tag_name,(obj.n_jets_max,),dtype=bool)
            # obj.buffer[self.tag_name] = np.full((obj.nevents,obj.n_jets_max),False,dtype=bool)
        return

    def _addFlagToBuffer(self,obj : 'JetFinder'):
        """
        Adds the ghost association tags to the buffer, for writing.
        """
        obj.buffer.set(self.tag_name,obj._i,[self.tags[i] for i in obj.jet_ordering])
        # embed_array_inplace([self.tags[i] for i in obj.jet_ordering],obj.buffer[self.tag_name][obj._i])

    def ModifyConstituents(self, obj : 'JetFinder'):
        return

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return
