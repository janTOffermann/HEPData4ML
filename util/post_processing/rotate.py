# The purpose of this code is to rotate four-vectors spatially, using two given
# truth-level 4-vectors as a reference (for how to rotate things).
# This may be useful for visualizations of the data.
import re
import subprocess as sub
import numpy as np
import h5py as h5
from util.post_processing.utils.rotations import DetermineRotation
import util.qol_utils.qol_util as qu

class Rotation:
    def __init__(self,truth_idx_1, truth_idx_2,rotation_keys=[".*Pmu.*"],exclusion_keys=[".*_cyl.*"],verbose=False):
        self.status = False
        self.SetVerbosity(verbose)
        self.SetTruthIndices(truth_idx_1,truth_idx_2)
        self.SetRotationKeyExpressions(rotation_keys)
        self.SetExclusionKeyExpressions(exclusion_keys)
        self.configurator = None
        self.print_prefix = '\n\tRotation'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tProducing rotated copies of four-vectors:'
        self.progress_bar_suffix = 'Complete'

    def RequiresIndexing(self):
        return False

    def SetVerbosity(self,flag):
        self.verbose = flag

    def SetH5EventFile(self,file):
        self.h5_file = file

    def SetTruthIndices(self,truth_idx_1,truth_idx_2):
        """
        Set the two truth indices used for the rotation -- the rotation will be done so that
        the first truth particle is centered, and the second is in a particular plane.
        """
        self.truth_idx_1 = truth_idx_1
        self.truth_idx_2 = truth_idx_2
        return

    def SetRotationKeyExpressions(self,expressions):
        """
        Sets a list of regular expressions used to determine which
        keys we will apply the rotations to.
        """
        self.rotation_key_expressions = expressions

    def SetExclusionKeyExpressions(self,expressions):
        """
        Sets a list of regular expressions used to determine which
        keys we will *not* apply rotations to -- only necessary if we
        want to exclude keys that will first be picked up by any of
        the patterns handed to SetRotationKeyExpressions().
        """
        self.exclusion_key_expressions = expressions

    def GetRotationMatrices(self):
        return self.rotation_matrices

    def GetRotatedData(self):
        return self.data_to_rotate # only will be rotated after running Process()

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def Initialize(self): # TODO: May want to consider chunking things and using a buffer? Memory usage will scale better for larger files.
        """
        Reads in the input HDF5 file, and places the required arrays in memory.
        """
        f = h5.File(self.h5_file,'r')
        self.keys = list(f.keys())
        self.truth_Pmu = f['truth_Pmu'][:]
        self.truth_Nobj = f['truth_Nobj'][:]
        self.nevents = self.truth_Pmu.shape[0]
        self.is_signal = f['is_signal'][:] # specifically to look for negative signal flags, which indicate events that are to be discarded (and will be lacking actual jets)

        # Now we will load the data corresponding to the data fields we are rotating.
        self.keys_to_rotate = []
        for key in self.keys:
            for expr in self.rotation_key_expressions:
                search = re.findall(expr,key)
                if(len(search) > 0):
                    if(key not in self.keys_to_rotate): self.keys_to_rotate.append(key)
                    break

        for expr in self.exclusion_key_expressions:
            searches = [re.findall(expr,key) for key in self.keys_to_rotate]
            mask = [len(x) == 0 for x in searches]
            self.keys_to_rotate = [b for a, b in zip(mask,self.keys_to_rotate) if a]

        self.keys_to_rotate.sort()
        self.data_to_rotate = {key:f[key][:] for key in self.keys_to_rotate}
        self.rotation_matrices = np.zeros((self.nevents,3,3),dtype=np.dtype('f8'))
        self.status = True
        f.close()

    def Process(self):
        self.Initialize()

        if(self.verbose):
            print('{}.Process: Number of events = {}'.format(self.print_prefix,self.nevents))
            qu.printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        for i in range(self.nevents):
            if(self.is_signal[i] < 0): continue
            v1 = self.truth_Pmu[i,self.truth_idx_1,:]
            v2 = self.truth_Pmu[i,self.truth_idx_2,:]
            R = DetermineRotation(v1,v2)
            self.rotation_matrices[i,:,:] = R # saving the rotation matrix itself

            # now apply the rotation matrix to the data
            for key in self.keys_to_rotate:
                self.data_to_rotate[key][i,...,1:4] = np.tensordot(self.data_to_rotate[key][i,...,1:4],R,axes=((-1,1))) # faster than einsum, and generalizes well!

            if(self.verbose):
                qu.printProgressBarColor(i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)
        return

    # Using a generic signature -- should consider making the various post-processors inherit from a single parent class!
    def __call__(self, delphes_file,h5_file,indices_file,output_file=None,verbose=None, copts=9, key=None):
        self.SetH5EventFile(h5_file)
        if(verbose is not None): self.SetVerbosity(verbose)
        self.Process()

        if(output_file is None): output_file = h5_file
        else: sub.check_call(['cp',h5_file,output_file])

        if(self.verbose):
            print('\tWriting rotation matrix and rotated four-vectors to {}.'.format(output_file))

        f = h5.File(output_file,'a')
        d = f.create_dataset('rotation_matrix',data=self.GetRotationMatrices(),compression='gzip',compression_opts=copts)
        for key,val in self.data_to_rotate.items():
            new_key = '{}_rot'.format(key)
            if(new_key in self.keys): continue
            d = f.create_dataset(new_key,data=val,compression='gzip',compression_opts=copts)
        f.close()
        return output_file