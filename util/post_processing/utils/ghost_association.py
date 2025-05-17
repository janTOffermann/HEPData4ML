import ROOT as rt
import h5py as h5
import numpy as np
import itertools

class GhostAssociator():

    def __init__(self,key,indices):


        self.vec_key = '{}.Pmu_cyl'.format(key)
        self.ghost_key = '{}.Ghost.Pmu'.format(key)
        self.indices = indices


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

        return

    def ModifyJets(self, obj):
        """
        This function will drop jets that aren't ghost-associated.
        """

        # Because of the indexing done in ModifyInitialization(), here
        # we know that the ghosts correspond to th first N particles.
        indices = np.arange(obj.input_collection_arrays[self.ghost_key].shape[1])

        # Fetch the jet constituents.
        obj._fetchJetConstituents()

        mask = np.full(len(obj.constituent_indices),False)

        for i,constituent_list in enumerate(obj.constituent_indices):
            for idx in indices:
                if(idx in constituent_list):
                    mask[i] = True
                    break

        obj.jets = list(itertools.compress(obj.jets,mask)) # NOTE: Not really necessary
        obj.jet_vectors = list(itertools.compress(obj.jet_vectors,mask))
        obj.jet_vectors_cyl = list(itertools.compress(obj.jet_vectors_cyl,mask))

        return

    def ModifyConstituents(self, obj):
        """
        Acts on obj.constituent_vectors etc.,
        to remove ghosts.
        """

        # Because of the indexing done in ModifyInitialization(), here
        # we know that the ghosts correspond to th first N particles.
        indices = np.arange(obj.input_collection_arrays[self.ghost_key].shape[1])

        for i in range(len(obj.constituent_vectors)):

            mask = [x not in indices for x in obj.constituent_indices[i]]
            obj.constituent_vectors[i] = list(itertools.compress(obj.constituent_vectors[i],mask))
            obj.constituent_vectors_cyl[i] = list(itertools.compress(obj.constituent_vectors_cyl[i],mask))
            obj.constituent_indices[i] = list(itertools.compress(obj.constituent_indices[i],mask))
        return