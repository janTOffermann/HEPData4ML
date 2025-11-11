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

class GhostAssociator(PostProcessorBase):
    """
    Ghost-associates jets with particles labeled by `key`,
    at the given particle indices (within each event).

    See: https://arxiv.org/abs/0802.1188 [JHEP 04 (2008) 005]
    """

    def __init__(self,key:str,indices,mode:str='filter',tag_name:Optional[str]=None):
        super().__init__()

        self.key = key
        self.vec_key_cyl = '{}.Pmu_cyl'.format(key)
        self.ghost_key = '{}.Ghost.Pmu'.format(key)
        self.indices = indices
        self.mode = mode
        assert self.mode in ['tag','filter']

        self.tag_name=tag_name
        self.tags = None
        self.name = 'GhostAssociator'
        self.print_prefix = '{}'.format(self.name)

        self.citations = {
            "GhostAssociation":
            """
@article{Cacciari:2008gn,
    author = "Cacciari, Matteo and Salam, Gavin P. and Soyez, Gregory",
    title = "{The Catchment Area of Jets}",
    eprint = "0802.1188",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "LPTHE-07-02",
    doi = "10.1088/1126-6708/2008/04/005",
    journal = "JHEP",
    volume = "04",
    pages = "005",
    year = "2008"
}
            """
        }

    def _makeGhosts(self,vecs:Union[list,np.ndarray], a:float=1.0e-10):
        """Make the ghosts, in both Cartesian and cylindrical."""
        result = (np.zeros(vecs.shape),np.zeros(vecs.shape))

        for i,vec in enumerate(vecs):
            v = rt.Math.PtEtaPhiMVector(a, vec[1], vec[2], a)
            result[0][i] = np.array([v.E(),v.Px(),v.Py(),v.Pz()])
            result[1][i] = np.array([v.Pt(), v.Eta(), v.Phi(), v.M()])
        return result

    def _input_check(self,obj:'JetFinder'):
        if(self.obj_name_input is not None):
            if(self.obj_name_input != obj.jet_name):
                print()
                self._print('Attempted to set input object name to: {}'.format(self.obj_name_input),level=1)
                self._print('This post-processor operates by modifying the base jet clustering itself,',level=1)
                self._print('so it doesn\'t take modified input. Will ignore this setting.',level=1)
                self._print('You might want to put this post-processor first in the chain;',level=1)
                self._print('that\'s effectively where it already operates.',level=1)

    def ModifyInitialization(self,obj:'JetFinder'):
        """
        This function will modify the initialization so that
        the ghost vectors are generated and loaded into memory.
        """
        self._input_check(obj)

        # NOTE: obj_name_input effectively unused, but want to pass info through
        if(self.obj_name_input is None):
            self.SetInputObjectName(obj.jet_name)
        self.obj_name_output = self.obj_name_input # doesn't modify jet -> pass through object name

        self.indices = np.atleast_1d(self.indices)

        # Fetch the 4-vector key, and make sure its data is loaded.
        # (If it's already loaded by the JetFinderi input_buffer,
        #  function won't do anything so no need to check).
        obj.input_buffer.read_branch(self.vec_key_cyl)

    def ModifyInputs(self,obj : 'JetFinder'):
        """
        Puts the ghost vectors corresponding with `key`
        on the bottom of obj.input_vecs.
        """

        ghost_input_vecs = np.array(obj.input_buffer[self.vec_key_cyl])[self.indices]
        ghost_vecs, ghost_vecs_cyl = self._makeGhosts(ghost_input_vecs)
        original_input_length = len(obj.input_vecs)

        obj.input_vecs    = np.vstack([obj.input_vecs,ghost_vecs])
        obj.input_vecs_cyl = np.vstack([obj.input_vecs_cyl,ghost_vecs_cyl]) # Needed for internal consistency throughout

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
            # obj._updateJetDictionary()
            # obj._jetsToVectors()
            # obj._fetchJetConstituents()
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
            self.tag_name = '{}.{}.GhostAssociated'.format(obj.jet_name,self.key) # NOTE: use of obj.jet_name
        if(self.tag_name not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.tag_name,ndim=1,dtype=bool)
        return

    def _addFlagToBuffer(self,obj : 'JetFinder'):
        """
        Adds the ghost association tags to the buffer, for writing.
        """
        obj.output_buffer.set(self.tag_name,obj._i,self.tags)
