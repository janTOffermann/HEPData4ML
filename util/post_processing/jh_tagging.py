# The purpose of this code is to apply the Johns Hopkins top tagger (arXiv:0806.0848 [hep-ph])
# to the jets in the dataset.
import subprocess as sub
import numpy as np
import h5py as h5
from util.post_processing.utils.jhtagger import JHTaggerSetup, JHTopTagger
from util.calcs import Calculator
import util.qol_utils.qol_util as qu

class JHTagger:
    def __init__(self,truth_idx=2, w_nobj_max = 120, verbose=False):
        self.status = False
        self.SetMaxWConstituents(w_nobj_max)
        self.SetTruthIndex(truth_idx)
        self.SetVerbosity(verbose)
        self.configurator = None
        self.calculator = None

        self.print_prefix = '\n\tJHTagger'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tRunning Johns Hopkins top tagger:'
        self.progress_bar_suffix = 'Complete'
        self.setup = None
        self.tagger = None

        self.error = False

    def RequiresIndexing(self):
        return False

    def SetVerbosity(self,flag):
        self.verbose = flag

    def SetMaxWConstituents(self,val):
        self.w_nobj_max = val

    def SetTruthIndex(self,val):
        self.w_index = val

    def SetH5EventFile(self,file):
        self.h5_file = file

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def Initialize(self): # TODO: May want to consider chunking things and using a buffer? Memory usage will scale better for larger files.
        """
        Reads in the input HDF5 file, and places the required arrays in memory.
        """

        # Set up the JHTagger.
        self.setup = JHTaggerSetup(self.configurator)
        self.setup.FullPreparation()
        self.tagger = JHTopTagger() # for now, just using default arguments

        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

        f = h5.File(self.h5_file,'r')
        self.Pmu = f['Pmu'][:]
        self.Nobj = f['Nobj'][:]
        self.nevents = self.Nobj.shape[0]
        self.is_signal = f['is_signal'][:] # specifically to look for negative signal flags, which indicate events that are to be discarded (and will be lacking actual jets)
        self.truth_Pmu = f['truth_Pmu'][:,self.w_index] # get the truth-level W -- used for calc'ing some of the things below

        self.jh_tag     = np.full(self.nevents,False,dtype=bool)
        self.jh_W_pmu   = np.zeros((self.nevents,4),dtype=np.dtype('f8'))
        self.jh_pt_pred = np.zeros(self.nevents,dtype=np.dtype('f8'))
        self.jh_m_pred  = np.zeros(self.nevents,dtype=np.dtype('f8'))
        self.jh_pt_res  = np.zeros(self.nevents,dtype=np.dtype('f8'))
        self.jh_m_res   = np.zeros(self.nevents,dtype=np.dtype('f8'))
        self.jh_W_pred_nobj = np.zeros(self.nevents,dtype=int)
        self.jh_W_pred_constituents = np.zeros((self.nevents,self.w_nobj_max,4),dtype=np.dtype('f8'))
        self.status = True
        f.close()

    def Process(self):
        self.Initialize()

        if(self.verbose):
            print('{}.Process: Number of events = {}'.format(self.print_prefix,self.nevents))
            qu.printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        for i in range(self.nevents):
            if(self.is_signal[i] < 0): continue
            self.tagger.TagEvent(self.Pmu[i])
            self.jh_tag[i] = self.tagger.GetStatus()
            error_flag = False

            if(self.jh_tag[i]):

                #TODO: Running into some errors on the UChicago Analysis Facility (CVMFS, with LCG_103).
                #      "cling JIT session error: Cannot allocate memory"
                #      Followed by crash in some line that uses a custom C++/ROOT library, like self.tagger funcs.
                #      As a workaround we can allow this to gracefully fail, and leave output empty.
                try: self.jh_W_pmu[i] = self.tagger.GetWCandidate()
                except:
                    self.jh_W_pmu[i,:] = np.nan # makes it obvious that something went wrong
                    error_flag = True

                try:
                    constituents = self.tagger.GetWCandidateConstituents()
                    self.jh_W_pred_nobj[i] = np.minimum(constituents.shape[0],self.w_nobj_max)
                    self.jh_W_pred_constituents[i,:self.jh_W_pred_nobj[i],:] = constituents
                except:
                    self.jh_W_pred_nobj[i] = -1
                    self.jh_W_pred_constituents[i,:,:] = np.nan
                    error_flag = True

                W_cyl      = self.calculator.EPxPyPzToPtEtaPhiM_single(*self.truth_Pmu[i])

                try:
                    W_pred_cyl = self.calculator.EPxPyPzToPtEtaPhiM_single(*self.jh_W_pmu[i])
                    self.jh_pt_pred[i] = W_pred_cyl[0]
                    self.jh_m_pred[i]  = W_pred_cyl[-1]
                    self.jh_pt_res[i] = W_pred_cyl[0]  / W_cyl[0]
                    self.jh_m_res[i]  = W_pred_cyl[-1] / W_cyl[-1]
                except:
                    W_pred_cyl = np.full(4,np.nan)
                    self.jh_pt_pred[i] = np.nan
                    self.jh_m_pred[i] = np.nan
                    self.jh_pt_res[i] = np.nan
                    self.jh_m_res[i] = np.nan
                    error_flag = True

            if(self.verbose):
                qu.printProgressBarColor(i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)
            if(error_flag and not self.error):
                self.error = True
                print(50*'=')
                print('\tWarning: Encountered error in JHTagger (e.g. a cling JIT memory allocation error).')
                print('\tWill continue, but fill JHTagger output with nan\'s.')
                print('\tThis will have to be remedied manually -- sorry for the inconvenience.')
                print(50*'=')
        return

    # Using a generic signature -- should consider making the various post-processors inherit from a single parent class!
    def __call__(self, delphes_file,h5_file,indices_file,output_file=None,verbose=None, copts=9, key=None):
        self.SetH5EventFile(h5_file)
        if(verbose is not None): self.SetVerbosity(verbose)
        self.Process()

        if(output_file is None): output_file = h5_file
        else: sub.check_call(['cp',h5_file,output_file])
        if(self.verbose):
            print('\tWriting JH tagging information to {}.'.format(output_file))
        f = h5.File(output_file,'a')
        d = f.create_dataset('jh_tag',              data=self.jh_tag,                compression='gzip',compression_opts=copts)
        d = f.create_dataset('jh_W_Pmu',            data=self.jh_W_pmu,              compression='gzip',compression_opts=copts)
        d = f.create_dataset('jh_W_constituent_Pmu',data=self.jh_W_pred_constituents,compression='gzip',compression_opts=copts)
        d = f.create_dataset('jh_W_Nobj',           data=self.jh_W_pred_nobj,        compression='gzip',compression_opts=copts)
        d = f.create_dataset('jh_m',                data=self.jh_m_pred,             compression='gzip',compression_opts=copts)
        d = f.create_dataset('jh_pt',               data=self.jh_pt_pred,            compression='gzip',compression_opts=copts)
        d = f.create_dataset('jh_m_resolution',     data=self.jh_m_res,              compression='gzip',compression_opts=copts)
        d = f.create_dataset('jh_pt_resolution',    data=self.jh_pt_res,             compression='gzip',compression_opts=copts)
        f.close()
        return output_file