import ROOT as rt
import numpy as np
import h5py as h5
import glob,sys,os,pathlib

from util.pileup.setup import PileupSetup

from util.qol_utils.pdg import DatabasePDG
from util.metadata.meta import MetaDataHandler
from util.config.config import Configurator
from util.hepmc.setup import HepMCSetup
from typing import Union, List, Tuple, Optional, TYPE_CHECKING

if(TYPE_CHECKING):
    setup = HepMCSetup(verbose=False)
    setup.PrepHepMC() # will download/install if necessary
    python_dir = setup.GetPythonDirectory()
    if(python_dir not in sys.path):
        sys.path = [setup.GetPythonDirectory()] + sys.path # prepend, to make sure we pick this one up first
    from pyHepMC3 import HepMC3 as hm

class PileupMixer:
    """
    This class performs on-the-fly mixing, to mix in events from some pileup HepMC3 file(s)
    into full "pileup events". These are placed in a new file, to be used as a "sidecar"
    with the HepMC3 files containing the main process.
    """

    def __init__(self, pileup_files:Optional[Union[str,list]]=None,rng_seed:int=1,mu_input:Optional[Union[str,List[Union[int,float]],Tuple[Union[int,float],Union[int,float]],np.ndarray]]=None, add_stable_only:bool=True, batch_size:int=100, contiguous_sampling:bool=True):

        # Immediately handle setup.
        # Out of an excess of caution,
        # make sure that HepMC3 is set up.

        self.setup_hepmc3 = HepMCSetup() # automatically runs

        self.setup = PileupSetup()
        self.setup.FullPreparation()

        self.mixer = rt.Pileup.PileupMixer()
        self.batch_size = batch_size
        self.mixer.SetBatchSize(self.batch_size)
        self.contiguous_sampling = contiguous_sampling
        self.mixer.SetUseContiguousSampling(self.contiguous_sampling)

        self.verbosity = 1
        self.name = 'PileupMixer'
        self.print_prefix = '{}'.format(self.name)

        self.require_pileup_input = True
        self.files = None
        self.SetPileupFiles(pileup_files)

        self.SetAddStableOnly(add_stable_only)

        self.selected_indices = None # transient storage for indices selected for a particular event
        self.allow_reuse = True
        self.pileup_indices = {} # accumulate the pileup event indices per event -> for output. Dictionary keys are filename
        self.number_non_pileup_particles = {} # accumulate number of non-pileup particles per event -- can be used later for indexing purposes
        self.n_particles_event = None # transient storage

        self.SetRNGSeed(rng_seed)
        self.beam_spot_sigma = None # will store the beam spot size
        self.SetBeamSpotSigma() # sets some sensible defaults
        self.do_phi_rotations = True
        self.phi_rotations_transient = []
        self.phi_rotations = {} # accumulate the phi rotations per event -> for output. Dictionary keys are filename

        # Set up the distribution for # of interactions per crossing
        self.available_mu_values = None
        self.mu_probabilities = None
        self.mu_input = mu_input
        self._init_mu_distribution()
        self.mu = None # transient storage for the current value of mu
        self.mu_values = {} # will store a list of all used mu values

        self.pileup_events = None # transient storage for pileup events

        self.event_buffer_size = 10 # how many combined events to store in memory, before flushing to output file
        #TODO: Need to fix buffer behavior -- the HepMC writers overwrite the output_file each time they're initialized, as opposed to appending

        self.indir = None
        self.outdir = None

        # for particle charge lookup
        self.pdg_database = DatabasePDG()

        self.metadata_handler = None
        self.configurator = None

        self.condor_flag = False
        self.condor_job_number = None
        self.condor_njobs = None
        self.condor_warning = False

    def SetAddStableOnly(self,flag:bool):
        self.filter_stable = flag
        self.mixer.SetStableOnly(self.filter_stable)

    def SetHTCondorInfo(self,flag:bool, job_number:int, njobs:int):
        # TODO: Deprecated/unused
        self.condor_flag = flag
        self.condor_job_number = job_number
        self.condor_njobs = njobs # the total number of jobs in this batch; can be useful for divvying up pileup events

    def SetMetadataHandler(self,handler:'MetaDataHandler'):
        self.metadata_handler = handler

    def SetConfigurator(self,configurator:'Configurator'):
        self.configurator = configurator

    def GetRNGSeed(self):
        return self.rng_seed

    def SetRNGSeed(self,rng_seed:int):
        if(rng_seed < 0):
            rng_seed = 0
        self.rng_seed = rng_seed
        self.rng = np.random.default_rng(self.rng_seed)

    def SetVerbosity(self,val:int):
        self.verbosity = val

    def SetInputDirectory(self,val:str):
        self.indir = val

    def SetOutputDirectory(self,val:str):
        self.outdir = val

    def SetAllowReuse(self, val:bool):
        self.allow_reuse = val

    def SetUsePhiRotations(self, val:bool):
        self.do_phi_rotations = val

    def SetPileupFiles(self,files:Union[str,list]=None):

        self.mixer.ClearPileupFiles()

        if(files is None):
            return
        if(isinstance(files,str)):
            self.files = glob.glob(files)
            if(len(self.files) == 0): # maybe this was not an absolute path, but a relative one (w.r.t. cwd)
                self.files = glob.glob('{}/{}'.format(os.getcwd(),files))
        else:
            self.files = files
            for i,file in enumerate(self.files):
                if(not pathlib.Path(file).exists()):
                    self.files[i] = '{}/{}'.format(os.getcwd(),file)
        self.files = sorted([str(pathlib.Path(x).absolute()) for x in self.files if pathlib.Path(x).exists()])

        for file in self.files:
            self.mixer.AddPileupFile(file)

        return

    def SetBeamSpotSigma(self,dt:float=0.16,dx:float=0.01,dy:float=0.01,dz:float=35.):
        """
        Sets the size of the beamspot (standard deviation) in (t,x,y,z), in nanoseconds or millimeters (as appropriate).
        The default spatial values are based on ATLAS Run 2: https://twiki.cern.ch/twiki/bin/view/AtlasPublic/BeamSpotPublicResults#Run_2_25ns_pp_Collisions_s_13_Te .
        The default time value is based on the Delphes ATLAS detector card (with pileup).
        """
        self.beam_spot_sigma = (dt,dx,dy,dz)
        self.mixer.SetBeamSpotSigma(*self.beam_spot_sigma)

    def _init_mu_from_file(self,file_name, hist_name=None):
        if(hist_name is None):
            hist_name = 'PileupOverlay_mu'
        try:
            f = rt.TFile(file_name,"READ")
            self.mu_distribution = f.Get(hist_name).Clone()
            self.mu_distribution.SetDirectory(0)
            f.Close()
        except:
            self._print("Error: failed to read histogram {} from file {}.".format(hist_name,file_name))
            self._print("Falling back on default mu distribution.")
            self.mu_input = None
            self._init_mu_distribution()

    def _init_mu_distribution(self):

        # Multiple kinds of information to parse.

        # None: use default distribution
        if(self.mu_input is None):
            # very approximate for Run 2,
            # see https://atlas.web.cern.ch/Atlas/GROUPS/DATAPREPARATION/PublicPlots/2018/DataSummary/figs/mu_2015_2018.png
            self.mu_input = [33.7,11.5]
            self._init_mu_distribution()
            return

        # String: interpret as a filepath, optionally with histogram name after a colon
        elif(isinstance(self.mu_input,str)):
            file_name = self.mu_input.split(':')[0]
            hist_name = None
            if(':' in self.mu_input):
                hist_name = self.mu_input.split(':')[-1]
            self._init_mu_from_file(file_name,hist_name)
            self.mixer.InitMuDistribution(self.mu_distribution)
            return

        elif(isinstance(self.mu_input,list) or isinstance(self.mu_input,tuple) or isinstance(self.mu_input,np.ndarray)):
            if(len(self.mu_input) == 2):
                self.mixer.InitMuDistribution(*self.mu_input) # Gaussian
                return

        self._print("mu_input not understood.")
        self._print("Falling back on default mu distribution.")
        self.mu_input = None
        self._init_mu_distribution()
        return

    def __call__(self,input_file:str,output_file:Optional[str]=None):
        from pyHepMC3 import HepMC3 as hm

        if('.hepmc.root' in input_file):
            input_file_extension = 'hepmc.root'
        else:
            input_file_extension = input_file.split('.')[-1]

        if(self.indir is not None):
            input_file = '{}/{}'.format(self.indir,input_file)

        if(output_file is None):
            output_file = input_file.replace('.{}'.format(input_file_extension),'.pileup.{}'.format(input_file_extension))
        elif(self.outdir is not None):
            output_file = '{}/{}'.format(self.outdir,output_file)

        self.mixer(input_file,output_file) # produces output_file, based on # of events in input_file
        return output_file

    def _record_pileup_info(self,output_file):

        #TODO: Rework this, need to do stuff on C++ side

        key = output_file.split('/')[-1] # remove any leading directory, just use filename (this is generally our convention)
        if(key not in self.pileup_indices.keys()):
            self.pileup_indices[key] = []
            self.phi_rotations[key] = []
            self.mu_values[key] = []
            self.number_non_pileup_particles[key] = []
        self.pileup_indices[key].append(self.selected_indices)
        self.mu_values[key].append(self.mu)
        self.number_non_pileup_particles[key].append(self.n_particles_event)

        # only record phi rotations if they were used -- no point in recording lots of zeros
        if(self.do_phi_rotations):
            self.phi_rotations[key].append(self.phi_rotations_transient)
        return

    def _migrate_pileup_info(self,old_key,new_key):
        """
        A bit of a messy function, a consequence of some file naming & I/O options.
        Simply renames a key in the info dictionaries.
        """
        for d in [self.mu_values,self.pileup_indices,self.phi_rotations,self.number_non_pileup_particles]:
            if(old_key in d.keys()):
                d[new_key] = d[old_key]
                del d[old_key]

    def Process(self,inputs,outputs=None):
       if(outputs is None):
        outputs = len(inputs) * [None]

        output_files = []
        for (input,output) in zip(inputs,outputs):
            output_file = self.__call__(input,output)
            output_files.append(output_file)

        #self.AddPileupInfoToHepMC3ROOT()

        self._writeMetadata()
        return output_files

    def _gaussian(self,x,mu,sig,A=None,require_positive=True):
        if(require_positive and x < 0.):
            return 0.
        if(A is None): A = 1. / (np.sqrt(2.0 * np.pi))
        return A * np.exp(-np.square((x - mu) / sig) / 2)

    def _writeMetadata(self):
        """
        Writes metadata to the metadata handler.
        """
        pileup_metadata = {}
        if(self.metadata_handler is None):
            self._print('Unable to write metadata; no handler was provided.')
            return

        # Here, the metadata from the pileup files themselves will be fetched.
        if(self.files is not None):
            for pileup_file in self.files:
                pileup_file_extension = pileup_file.split('.')[-1]
                if(pileup_file_extension.lower() != 'root'):
                    self._print('Warning: Cannot read in metadata from pileup file {} .'.format(pileup_file))
                    continue
                key = pileup_file.split('/')[-1]
                pileup_metadata[key] = self.metadata_handler.ReadMetaDataFromROOTFile(pileup_file)
            self.metadata_handler.AddElement('Metadata.Pileup.InputMetadata',pileup_metadata)

        # Also add info on the mu distribution that was used
        self.mu_distribution = self.mixer.GetMuDistribution()
        mu_bins = np.array([self.mu_distribution.GetBinLowEdge(x+1) for x in range(self.mu_distribution.GetNbinsX())]) # left edges)
        mu_weights = np.array([self.mu_distribution.GetBinContent(x+1) for x in range(self.mu_distribution.GetNbinsX())])
        self.metadata_handler.AddElement('Metadata.Pileup.MuDistribution.BinEdgesLeft',mu_bins)
        self.metadata_handler.AddElement('Metadata.Pileup.MuDistribution.BinContents',mu_weights)
        return

    def GetMuValues(self):
        return self.mu_values

    def GetPileupInfo(self):
        data = {}
        data['Pileup.Mu'] = self.mu_values
        data['Pileup.Index'] = self.pileup_indices # index w.r.t. input pileup collection
        if(self.do_phi_rotations):
            data['Pileup.PhiRotation'] = self.phi_rotations
        return data

    def AddPileupInfoToHepMC3ROOT(self):
        """
        If producing output HepMC3/ROOT files, we actually add the
        pileup information (mu, indices etc.) to the HepMC3/ROOT file itself,
        in a separate tree. Note that we currently add this to the HDF5 files
        in an independent manner, but this may be helpful for debugging purposes
        and in the future could be the way to propagate the information to the
        final n-tuple.
        """

        # NOTE: With our current structure of GetPileupInfo() output, it would be
        # most natural to iterate over branches, and iterate over files in a nested loop.
        # However, for I/O reasons it's nicer to iterate over the files. Given how the
        # output is structured, it is safe to do this as we can assume each nested dict has
        # the same keys -- but this is generally kind of clunky code! -Jan
        data = self.GetPileupInfo()
        files = data[list(data.keys())[0]].keys() # NOTE: These are filenames *without* directories (as usual in this cocde).

        for i,file in enumerate(files):

            filename_full = '{}/{}'.format(self.outdir,file)

            f = rt.TFile(filename_full,'UPDATE')
            t = rt.TTree('PileupInfo','Pileup information from {}'.format(self.name))

            buffer_dict = {}

            for branch in data.keys():
                value = data[branch][file]

                is_nested_list = False
                if(isinstance(value[0],list)):
                    is_nested_list = True
                elif(isinstance(value[0],np.ndarray)):
                    is_nested_list = True

                if(is_nested_list):
                    #NOTE: Assuming dim=2

                    if(isinstance(value[0][0],int)): # NOTE: Would break if length==0. Would this ever happen?
                        buffer_dict[branch] = rt.std.vector('int')()
                    else:
                        buffer_dict[branch] = rt.std.vector('double')()
                    t.Branch(branch,buffer_dict[branch])

                else:
                    # simple 1D array
                    if(isinstance(value[0],int)):
                        buffer_dict[branch] = np.zeros(1,dtype=int)
                        t.Branch(branch,buffer_dict[branch],'{}/I'.format(branch))

                    else:
                        buffer_dict[branch] = np.zeros(1,dtype=float)
                        t.Branch(branch,buffer_dict[branch],'{}/D'.format(branch))

            # branch buffers are created, time to fill
            nentries = len(data[list(data.keys())[0]][file])
            for j in range(nentries):
                for branch in buffer_dict.keys():

                    if(isinstance(buffer_dict[branch], np.ndarray)):
                        buffer_dict[branch][0] = data[branch][file][j]
                    else:
                        for k,entry in enumerate(data[branch][file][j]):
                            buffer_dict[branch].push_back(entry)

                t.Fill()
            t.Write()
            f.Close()
        return

    def AddPileupInfoToH5(self,h5_file,file_key, cwd=None,copts=9):
        if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)

        f = h5.File(h5_file,'r+')
        keys = list(f.keys())
        # nevents = f[keys[0]].shape[0]
        data = self.GetPileupInfo()

        for key,value_dict in data.items():
            # for lists of (variable-length) arrays, need to embed these in some fixed length array

            value = value_dict[file_key]

            is_nested_list = False
            if(isinstance(value[0],list)):
                is_nested_list = True
            elif(isinstance(value[0],np.ndarray)):
                is_nested_list = True

            if(is_nested_list):
                # create an array
                max_length = np.max([len(x) for x in value])
                value_array = np.array([np.pad(row, (0, max_length-len(row))) for row in value])
                f.create_dataset(key,data=value_array,compression='gzip',compression_opts=copts)
            else:
                f.create_dataset(key,data=value,compression='gzip',compression_opts=copts)

        # With pileup information added, we can actually determine which entries in StableTruthParticles
        # came from pileup, so we can now do some post-processing to add that in
        particle_index_keys = [x for x in f.keys() if 'HepMC3Index' in x]

        for key in particle_index_keys:
            key_prefix = '.'.join(key.split('.')[:-1])
            pileup_flag_key = '{}.IsPileup'.format(key_prefix)
            is_pileup = np.array(f[key][:] > np.array(self.number_non_pileup_particles[file_key])[:,np.newaxis],dtype=bool) # indices use 1-indexing
            f.create_dataset(pileup_flag_key,data=is_pileup,dtype=bool,compression='gzip',compression_opts=copts)

        f.close()

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return