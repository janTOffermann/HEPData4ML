import sys,os,importlib,pathlib
import pathlib
import numpy as np
# from config.config import config

def GetConfigDictionary(config_file):
    config_name = config_file.split('/')[-1].split('.')[0]
    spec = importlib.util.spec_from_file_location("config.{}".format(config_name),config_file)
    config = importlib.util.module_from_spec(spec)
    sys.modules['config.{}'.format(config_name)] = config
    spec.loader.exec_module(config)
    return config.config

def Bool2String(bool):
    if(bool): return 'on'
    return 'off'

# Function for fetching the config.py file as a string.
def GetConfigFileContent(config_file):
    with open(config_file,'r') as f:
        lines = f.readlines()
    lines = [x.strip().strip('\n') for x in lines]
    return lines

class Configurator:
    def __init__(self,config_dictionary):
        self.generation=config_dictionary['generation']
        self.simulation=config_dictionary['simulation']
        self.reconstruction=config_dictionary['reconstruction']

        self.config = {}
        self.config['generation'] = self.generation
        self.config['simulation'] = self.simulation
        self.config['reconstruction'] = self.reconstruction

        # Will also use this for some print statement bookkeeping, since it is conveniently passed around.
        self.print_fastjet = True
        self.print_delphes = True

        self.bad_filepath_keys = []
        self._check_paths()

    def _check_paths(self):
        """
        Checks particular entries, which are filepaths.
        These could be problematic if they are given as relative paths,
        especially if using a running option where the code isn't copied
        inside the running directory.
        """
        self.status = True
        path_keys = ['delphes_dir','delphes_card','fastjet_dir']

        path_keys = {
            'simulation':['delphes_dir','delphes_card'],
            'reconstruction':['fastjet_dir']
        }

        for key,key_list in path_keys.items():
            for key2 in key_list:
                status = True
                if(key2 in self.config[key].keys()):
                    val = self.config[key][key2]
                    if(val is None):
                        continue # for some keys, it is OK to be "None"
                    else:
                        status = pathlib.Path(val).exists()
                if(not status):
                    self.status = False
                    print('Warning: Configuration entry [{}] = "{}", which doesn\'t appear to be a valid filepath.'.format(key2,val))
                    self.bad_filepath_keys.append('{}:{}'.format(key,key2))
                    return

    def CorrectFilepaths(self,dir):
        """
        Attempt to correct filepaths.
        """
        for key_pair in self.bad_filepath_keys:
            key,key2 = key_pair.split(':')
            old_filepath = self.config[key][key2]
            new_filepath = '{}/{}'.format(dir,old_filepath)
            print('\tAttempting to correct: {} -> {}'.format(old_filepath,new_filepath))
            self.config[key][key2] = new_filepath
        self.bad_filepath_keys = []
        self._check_paths()
        return

    def GetStatus(self):
        return self.status

    def GetHepMCFormat(self):
        return self.config['generation']['hepmc_format']

    def GetPythiaConfigFile(self,filename=None):
        """
        Fetch the Pythia configuration file.
        """
        if(filename is None):
            proc = self.config['generation']['proc']
            if('.' not in proc.split('/')[-1]):
                proc += '.txt'
            try:
                template_file = proc
                assert pathlib.Path(template_file).exists()
            except:
                path_to_this = os.path.dirname(os.path.realpath(__file__))
                template_file = '{}/pythia_templates/{}'.format(path_to_this,proc)
        else: template_file = filename
        if(not pathlib.Path(template_file).exists()):
            print('Error: Template file {} not found.'.format(template_file))
            assert(False) # TODO: Nicer error handling
        return template_file

    def GetPythiaConfigFileContents(self,filename=None):
        filepath = self.GetPythiaConfigFile(filename)
        with open(filepath,'r') as f:
            contents = f.readlines()
        return ''.join(contents)

    def GetPythiaVerbosity(self):
        return self.config['generation']['verbose']

    # Get the Pythia configuration, for a single pT bin, as a dictionary.
    # This dictionary will only control a few settings (MPI, ISR/FSR etc.),
    # while the process settings from the 'proc' entry of config will be
    # passed separately.
    def GetPythiaConfig(self,pt_min, pt_max,verbose=False):
        pythia_config = {}
        # Tack on a few more things to the pythia configuration
        pythia_config['HadronLevel:all'] = Bool2String(self.config['generation']['hadronization'])
        pythia_config['PartonLevel:MPI'] = Bool2String(self.config['generation']['mpi'])
        pythia_config['PartonLevel:ISR'] = Bool2String(self.config['generation']['isr'])
        pythia_config['PartonLevel:FSR'] = Bool2String(self.config['generation']['fsr'])
        pythia_config['Random:setSeed'] = 'on'
        pythia_config['Random:seed'] = self.config['generation']['rng']

        # Add the phase space stuff.
        pythia_config['PhaseSpace:pTHatMin'] = pt_min
        pythia_config['PhaseSpace:pTHatMax'] = pt_max

        if(not verbose):
            pythia_config['Print:quiet'] = 'on' # avoid printing reams of info
            pythia_config['Stat:showProcessLevel'] = 'off'
            pythia_config['Stat:showErrors'] = 'off'
        else:
            pythia_config['Print:quiet'] = 'off'
            pythia_config['Stat:showProcessLevel'] = 'on'
            pythia_config['Stat:showErrors'] = 'on'
            pythia_config['Next:numberShowProcess'] = '1'
            pythia_config['Next:numberShowEvent'] = '1'
        return pythia_config

    def GetPythiaRNGSeed(self):
        return self.config['generation']['rng']

    def GetParticleSelection(self):
        return self.config['reconstruction']['particle_selection']

    # def GetJetConfig(self):
    #     return_dict = {}
    #     for key in ['jet_radius','jet_min_pt','jet_max_eta','jet_n_par','jet_selection']: return_dict[key] = self.config['reconstruction'][key]
    #     return return_dict

    def GetSimulationType(self):
        return self.config['simulation']['type'].lower()

    def SetDelphesConfig(self,value=True):
        """
        Allows us to (re)set the Delphes flag.
        """
        if(type(value) != bool):
            print('Configurator.SetDelphesConfig(): Input {} not understood, setting delphes=False.'.format(value))
            value = False
        self.config['simulation']['type'] = value # delphes config
        return

    def GetDelphesCard(self):
        return self.config['simulation']['delphes_card']

    def GetNPars(self):
        # some fanciness for the number of delphes objects, since it can be a list
        if(len(np.asarray(self.config['reconstruction']['n_delphes'])) == 1):
            self.config['reconstruction']['n_delphes'] = np.full(len(self.config['simulation']['delphes_output']), np.asarray(self.config['reconstruction']['n_delphes'])[0],dtype=int)
        return_dict = {}
        for key in ['n_truth','n_stable','n_delphes']: return_dict[key] = self.config['reconstruction'][key]
        return return_dict

    def GetSignalFlag(self):
        return self.config['reconstruction']['signal_flag']

    def GetSplitSeed(self):
        return self.config['reconstruction']['split_seed']

    def GetPostProcessing(self):
        return self.config['reconstruction']['post_processing']

    def GetEventFilter(self):
        return self.config['reconstruction']['event_filter']

    def GetEventFilterFlag(self):
        return self.config['reconstruction']['event_filter_flag']

    def GetDelphesDirectory(self):
        return self.config['simulation']['delphes_dir']

    def SetDelphesDirectory(self,val):
        self.config['simulation']['delphes_dir'] = val

    def GetFastjetDirectory(self):
        return self.config['reconstruction']['fastjet_dir']

    def GetPrintFastjet(self):
        return self.print_fastjet

    def GetPrintDelphes(self):
        return self.print_delphes

    def SetPrintFastjet(self,val):
        self.print_fastjet = val

    def SetPrintDelphes(self,val):
        self.print_delphes = val

    # def GetUseVectorCalcs(self): # temporarily (?) disabling
    #     return self.config['use_vectorcalcs']

    def GetDelphesObjects(self):
        key = 'delphes_output'
        if(key in self.config['simulation'].keys()):
            return self.config['simulation'][key]
        else: return ['Tower'] # some default result, maybe good for backwards compatibility!