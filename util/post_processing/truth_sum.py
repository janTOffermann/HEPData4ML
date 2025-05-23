import subprocess as sub
import numpy as np
import h5py as h5
from util.calcs import Calculator
import util.qol_utils.qol_util as qu
from util.config import Configurator

class TruthParticleSumNearJet:
    """
    The purpose of this code is to get the four-vector sum of a set of particles
    within some distance of the jet centroid. Specifically, one gives a set of
    truth particle indices (w.r.t. whatever truth particles have been saved in truth_Pmu),
    and this code will sum whichever ones are within some specified dR of jet_Pmu.
    For example, one might use this to find the four-vector sum of the daughter particles
    of a W-boson from top quark decay, that fall within the top quark jet (the daughters will
    have to be saved in the truth_Pmu array).
    """

    def __init__(self,truth_indices=[],jet_distance=0.8,verbose=False):
        self.status = False
        self.SetVerbosity(verbose)
        self.SetTruthIndices(truth_indices)
        self.SetJetDistance(jet_distance)
        self.configurator = None
        self.calculator = None
        self.calculator = Calculator(use_vectorcalcs=False) # for now, default to false instead of having to take in the configurator

        self.print_prefix = '\n\tTruthParticleSumNearJet'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tSumming selected truth particles within ∆R < {:.1f} of jet:'.format(self.jet_distance)
        self.progress_bar_suffix = 'Complete'
        self.pmu_sum = None

    def RequiresIndexing(self):
        return False

    def SetVerbosity(self,flag):
        self.verbose = flag

    def SetH5EventFile(self,file):
        self.h5_file = file

    def SetTruthIndices(self,truth_indices):
        self.truth_indices = np.sort(np.array(truth_indices,dtype=int)) # need to sort for h5py I/O

    def SetJetDistance(self,val):
        """
        Sets the maximum distance, in (eta,phi), between particles we'll
        sum and the jet centroid. Any truth particles picked up by the supplied
        truth indices, but outside of this distance, will be omitted from the sum.
        """
        self.jet_distance = val
        self.progress_bar_prefix = '\tSumming selected truth particles within dR = {:.1f} of jet:'.format(self.jet_distance)

    def GetPmuSum(self):
        return self.pmu_sum

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def Initialize(self): # TODO: May want to consider chunking things and using a buffer? Memory usage will scale better for larger files.
        """
        Reads in the input HDF5 file, and places the required arrays in memory.
        """
        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())
        f = h5.File(self.h5_file,'r')
        self.truth_Pmu = f['truth_Pmu'][:,self.truth_indices,:] # can have some zeros, if a particular event doesn't have a particle at an index in self.truth_indices
        self.truth_Nobj = f['truth_Nobj'][:]
        self.jet_Pmu_cyl = f['jet_Pmu_cyl'][:]
        self.nevents = self.jet_Pmu_cyl.shape[0]
        self.is_signal = f['is_signal'][:] # specifically to look for negative signal flags, which indicate events that are to be discarded (and will be lacking actual jets)
        self.pmu_sum = np.zeros((self.nevents,4),dtype=np.dtype('f8'))
        self.status = True
        f.close()

    def Process(self):
        self.Initialize()

        if(self.verbose):
            print('{}.Process: Number of events = {}'.format(self.print_prefix,self.nevents))
            qu.printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        for i in range(self.nevents):
            if(self.is_signal[i] < 0): continue
            truth_Pmu = self.truth_Pmu[i][self.truth_indices < self.truth_Nobj[i]] # masking to get rid of any empty entries
            jet_Pmu = self.jet_Pmu_cyl[i]
            truth_Pmu_cyl = self.calculator.PxPyPzEToPtEtaPhiM(*np.roll(truth_Pmu,-1,axis=1).T)
            dr2 = self.calculator.DeltaR2Vectorized(truth_Pmu_cyl[:,1:3],np.expand_dims(jet_Pmu[1:3],0)).flatten()
            mask = dr2 <= np.square(self.jet_distance)
            self.pmu_sum[i,:] = np.sum(truth_Pmu[mask,:],axis=0)
            if(self.verbose):
                qu.printProgressBarColor(i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)
        return

    # Using a generic signature -- should consider making the various post-processors inherit from a single parent class!
    def __call__(self, delphes_file,h5_file,indices_file,output_file=None,verbose=None, copts=9, key='contained_daughter_sum_Pmu'):
        self.SetH5EventFile(h5_file)
        if(verbose is not None): self.SetVerbosity(verbose)
        self.Process()

        if(output_file is None): output_file = h5_file
        else: sub.check_call(['cp',h5_file,output_file])

        if(self.verbose):
            print('\tWriting {} to {}.'.format(key,output_file))
        f = h5.File(output_file,'a')
        d = f.create_dataset(key,data=self.GetPmuSum(),compression='gzip',compression_opts=copts)
        f.close()
        return output_file

class TruthParticleSumByPdgID:
    """
    Produce a four-momentum sum of all truth-particles in an event, that have a PDG ID
    contained in the provided list. Produces output in both Cartesian and cylindrical.
    """

    def __init__(self,truth_pdg=[],compute_cartesian=True,compute_cylindrical=True, compute_mass=False, verbose=False):
        self.status = False
        self.SetVerbosity(verbose)
        self.SetTruthPdg(truth_pdg)
        self.configurator = None
        self.calculator = None
        self.calculator = Calculator(use_vectorcalcs=False) # for now, default to false instead of having to take in the configurator

        self.print_prefix = '\n\tTruthParticleSumByPdgID'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tSumming selected truth particles:'
        self.progress_bar_suffix = 'Complete'
        self.pmu_sum = None
        self.pmu_sum_cyl = None
        self.mass = None

        self.compute_cartesian = compute_cartesian
        self.compute_cylindrical = compute_cylindrical
        self.compute_mass = compute_mass

    def RequiresIndexing(self):
        return False

    def SetVerbosity(self,flag):
        self.verbose = flag

    def SetH5EventFile(self,file):
        self.h5_file = file

    def SetTruthPdg(self,truth_pdg):
        self.truth_pdg_list = np.array(truth_pdg,dtype=int)

    def GetPmuSum(self):
        return self.pmu_sum

    def GetPmuSumCyl(self):
        return self.pmu_sum_cyl

    def GetMass(self):
        return self.mass

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def Initialize(self): # TODO: May want to consider chunking things and using a buffer? Memory usage will scale better for larger files.
        """
        Reads in the input HDF5 file, and places the required arrays in memory.
        """
        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())
        f = h5.File(self.h5_file,'r')
        self.truth_Pmu = f['truth_Pmu'][:] # can have some zeros, if a particular event doesn't have a particle at an index in self.truth_indices
        self.truth_Pdg = f['truth_Pdg'][:] # can have some zeros, if a particular event doesn't have a particle at an index in self.truth_indices

        self.truth_Nobj = f['truth_Nobj'][:]
        self.nevents = self.truth_Pmu.shape[0]
        self.is_signal = f['is_signal'][:] # specifically to look for negative signal flags, which indicate events that are to be discarded
        self.pmu_sum = np.zeros((self.nevents,4),dtype=np.dtype('f8'))
        self.pmu_sum_cyl = np.zeros((self.nevents,4),dtype=np.dtype('f8'))
        self.mass = np.zeros(self.nevents,dtype=np.dtype('f8'))

        self.status = True
        f.close()

    def Process(self):
        self.Initialize()

        if(self.verbose):
            print('{}.Process: Number of events = {}'.format(self.print_prefix,self.nevents))
            qu.printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        for i in range(self.nevents):
            if(self.is_signal[i] < 0): continue
            mask = [x in self.truth_pdg_list for x in self.truth_Pdg[i]]
            truth_Pmu = self.truth_Pmu[i][mask]
            self.pmu_sum[i,:] = np.sum(truth_Pmu,axis=0)

            if(self.compute_cylindrical):
                self.pmu_sum_cyl[i,:] = self.calculator.EPxPyPzToPtEtaPhiM_single(*self.pmu_sum[i,:])

            if(self.compute_mass):
                self.mass[i] = np.sqrt( np.square(self.pmu_sum[i,0]) - np.dot(self.pmu_sum[i,1:],self.pmu_sum[i,1:]) ) # TODO: should probably have formula in calculator

            if(self.verbose):
                qu.printProgressBarColor(i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)
        return

    # Using a generic signature -- should consider making the various post-processors inherit from a single parent class!
    def __call__(self, delphes_file,h5_file,indices_file,output_file=None,verbose=None, copts=9, key='truth_Pmu_sum'):
        self.SetH5EventFile(h5_file)
        if(verbose is not None): self.SetVerbosity(verbose)
        self.Process()

        if(output_file is None): output_file = h5_file
        else: sub.check_call(['cp',h5_file,output_file])

        key_cyl = key + '_cyl'
        key_mass = key + '_mass'

        f = h5.File(output_file,'a')
        if(self.compute_cartesian):
            if(self.verbose):
                print('\tWriting {} to {}.'.format(key,output_file))
            d = f.create_dataset(key,data=self.GetPmuSum(),compression='gzip',compression_opts=copts)

        if(self.compute_cylindrical):
            if(self.verbose):
                print('\tWriting {} to {}.'.format(key_cyl,output_file))
            d = f.create_dataset(key_cyl,data=self.GetPmuSumCyl(),compression='gzip',compression_opts=copts)

        if(self.compute_mass):
            if(self.verbose):
                print('\tWriting {} to {}.'.format(key_mass,output_file))
            d = f.create_dataset(key_mass,data=self.GetMass(),compression='gzip',compression_opts=copts)
        f.close()
        return output_file

