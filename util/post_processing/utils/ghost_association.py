import itertools
import ROOT as rt
import h5py as h5
import numpy as np
from util.calcs import embed_array_inplace

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

    def ModifyInitialization(self,obj):
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

    def ModifyInputs(self,obj):
        """
        Puts the ghost vectors corresponding with `key`
        on the top of obj.input_vecs.
        """
        ghost_vecs = obj.input_collection_arrays[self.ghost_key][obj._i] # ith event
        obj.input_vecs = np.vstack([ghost_vecs,obj.input_vecs])

        for i in range(len(obj.input_vecs)):
            ghost_dict = {'GhostAssociation:Ghost':(i < len(ghost_vecs))}
            obj.AddUserInfo(i,ghost_dict)
        return

    def ModifyJets(self, obj):
        """
        This function will do one of two things, depending on self.mode:
        - self.mode == 'filter': Drop jets that aren't ghost-associated.
        - self.mode =='tag': Write outputs indicating whether or not the jet is ghost-associated.
        This function will also make sure that ghosts in the jets' constituent lists are removed.
        """
        import fastjet as fj # NOTE: In practice, fastjet will have been initialized already by JetFinder. Can similarly do this in Softdrop

        # Fetch the jet constituents. This fills obj.constituent_indices
        obj._fetchJetConstituents()

        self.tags = np.full(len(obj.jets),False)

        for i,jet in enumerate(obj.jets):
            ghost_mask = np.array([x.python_info()['GhostAssociation:Ghost'] for x in jet.constituents()])
            self.tags[i] = np.sum(ghost_mask) > 0

            # For the jets with ghosts, modify them to remove the ghost --
            # we don't want to pass it to anyfurther steps.
            if(self.tags[i]):
                obj.jets[i] = fj.join([x for x in list(itertools.compress(list(jet.constituents()),~ghost_mask))])

        if(self.mode=='filter'):
            obj.jets            = list(itertools.compress(obj.jets,self.tags           ))
            obj._jetsToVectors() # refresh the jet vectors, to deal with the filtering
            obj._fetchJetConstituents()
        # else:
        #     for i,entry in enumerate(self.tags):
        #         self._addFlagToBuffer(obj,i,entry)

        return

    def ModifyWrite(self,obj):
        if(self.mode=='filter'):
            return # do nothing
        else:
            self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
            self._addFlagToBuffer(obj)

    def _initializeBuffer(self,obj):
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
            obj.buffer[self.tag_name] = np.full((obj.nevents,obj.n_jets_max),False,dtype=bool)
        return

    def _addFlagToBuffer(self,obj):
        """
        Adds the JH tags to the buffer, for writing.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        #TODO: Is there a cleaner way to handle the pT sorting?
        #      Is this current implementation robust?
        embed_array_inplace(self.tags[obj.pt_sorting],obj.buffer[self.tag_name][obj._i])

    def ModifyConstituents(self, obj):
        return

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return
