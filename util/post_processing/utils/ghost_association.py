import sys, itertools
import ROOT as rt
import h5py as h5
import numpy as np
# from util.fastjet import JetFinderBase

# class GhostAssociatorv2():
#     """
#     Ghost-associates jets with particles labeled by `key`,
#     at the given particle indices (within each event).

#     See: https://arxiv.org/abs/0802.1188 [JHEP 04 (2008) 005]
#     """

#     def __init__(self,key,indices):

#         self.vec_key = '{}.Pmu_cyl'.format(key)
#         self.ghost_key = '{}.Ghost.Pmu'.format(key)
#         self.indices = indices
#         self.input_collection_arrays = {}

#         self.jet_finder = JetFinderBase()

#     def MakeGhosts(self,vecs, a=1.0e-10):

#         result = np.zeros(vecs.shape)
#         for i,vec in enumerate(vecs):
#             v = rt.Math.PtEtaPhiMVector(a, vec[1], vec[2], a)
#             result[i] = np.array([v.E(),v.Px(),v.Py(),v.Pz()])
#         return result

#     def ModifyInitialization(self,obj):
#         """
#         This function will modify the initialization so that
#         the ghost vectors are generated and loaded into memory.
#         """

#         self.jet_finder.SetConfigurator(obj.configurator)
#         self.jet_finder.SetRadius(obj.radius)
#         self.jet_finder.jet_algorithm_name = obj.jet_algorithm_name

#         # Initialize fastjet.
#         self.jet_finder._initialize_fastjet()

#         self.indices = np.atleast_1d(self.indices)
#         # Fetch the 4-vector key, and make sure its data is loaded.
#         # Note the use of cylindrical coordinates!
#         vecs = None
#         if(self.vec_key not in obj.input_collection_arrays.keys()):
#             # Read the necessary keys in.
#             f = h5.File(obj.h5_file,'r')

#             # Cylindrical coordinates four-momenta
#             vecs = f[self.vec_key][:][:,self.indices] # only loads the data needed -- indexing already done here

#             f.close()
#         else:
#             vecs = obj.input_collection_arrays[self.vec_key]

#         # Now convert these to ghosts: send pT and m -> 0.
#         # These are in Cartesian (E,px,py,pz)
#         ghost_vecs = np.array([self.MakeGhosts(v) for v in vecs])
#         self.input_collection_arrays[self.ghost_key] = ghost_vecs

#     def ModifyInputs(self,obj):
#         # NOTE: Doesn't do anything, we don't want to inject the ghosts into the
#         # JetFinder inputs because these could cause some interplay with other processors
#         # and lead to unintended effects. (Namely, we'll have a hard time removing the ghosts
#         # later, if their index in the jet is changed by another processor).
#         return

#     def ModifyJets(self, obj):
#         """
#         This function will re-run jet clustering with ghosts,
#         to determine which jets are ghost-associated, and then
#         filter the original jets to only keep those that are.
#         """

#         # Because of the indexing done in ModifyInitialization(), here
#         # we know that the ghosts correspond to th first N particles.
#         indices = np.arange(obj.input_collection_arrays[self.ghost_key].shape[1])

#         # Fetch the jet constituents.
#         obj._fetchJetConstituents()

#         mask = np.full(len(obj.constituent_indices),False)

#         for i,constituent_list in enumerate(obj.constituent_indices):
#             for idx in indices:
#                 if(idx in constituent_list):
#                     mask[i] = True
#                     break

#         obj.jets            = list(itertools.compress(obj.jets,mask           ))
#         obj.jet_vectors     = list(itertools.compress(obj.jet_vectors,mask    ))
#         obj.jet_vectors_cyl = list(itertools.compress(obj.jet_vectors_cyl,mask))

#         return


class GhostAssociator():
    """
    Ghost-associates jets with particles labeled by `key`,
    at the given particle indices (within each event).

    See: https://arxiv.org/abs/0802.1188 [JHEP 04 (2008) 005]
    """

    def __init__(self,key,indices):

        self.vec_key = '{}.Pmu_cyl'.format(key)
        self.ghost_key = '{}.Ghost.Pmu'.format(key)
        self.indices = indices

    # def _initialize_fastjet(self):
    #     """
    #     Based on JetFinderBase._initialize_fastjet().
    #     Not sure if this is really needed in practice.
    #     """
    #     self.jet_finder._setupFastJet()
    #     sys.path.append(self.jet_finder.fastjet_dir)
    #     import fastjet as fj


    def MakeGhosts(self,vecs, a=1.0e-10):

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
        ghost_vecs = np.array([self.MakeGhosts(v) for v in obj.input_collection_arrays[self.vec_key]])

        obj.input_collection_arrays[self.ghost_key] = ghost_vecs

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
        This function will drop jets that aren't ghost-associated.
        """
        import fastjet as fj

        # Fetch the jet constituents. This fills obj.constituent_indices
        obj._fetchJetConstituents()

        mask = np.full(len(obj.jets),False)

        for i,jet in enumerate(obj.jets):
            ghost_mask = np.array([x.python_info()['GhostAssociation:Ghost'] for x in jet.constituents()])
            mask[i] = np.sum(ghost_mask) > 0

            # For the jets we keep, modify them to remove the ghost -- we don't want to pass it to any
            # further steps.
            if(mask[i]):
                obj.jets[i] = fj.join([x for x in list(itertools.compress(list(jet.constituents()),~ghost_mask))])

        obj.jets            = list(itertools.compress(obj.jets,mask           ))
        obj.jet_vectors     = list(itertools.compress(obj.jet_vectors,mask    ))
        obj.jet_vectors_cyl = list(itertools.compress(obj.jet_vectors_cyl,mask))
        return


    def ModifyConstituents(self, obj):
        return
