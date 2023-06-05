import sys,os
import pathlib
# from config.config import config

def Bool2String(bool):
    if(bool): return 'on'
    return 'off'

# Function for fetching the config.py file as a string.
def GetConfigFileContent():
    config_file = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/../config/config.py') # assuming relative location is fixed
    with open(config_file,'r') as f:
        lines = f.readlines()
    lines = [x.strip().strip('\n') for x in lines]
    return lines

class Configurator:
    def __init__(self,config_dictionary):
        self.config=config_dictionary

        # Will also use this for some print statement bookkeeping, since it is conveniently passed around.
        self.print_fastjet = True
        self.print_delphes = True

    def GetPythiaConfigFile(self,filename=None):
        if(filename is None):
            path_to_this = os.path.dirname(os.path.realpath(__file__))
            proc = self.config['proc']
            template_file = '{}/pythia_templates/{}.txt'.format(path_to_this,proc)
        else: template_file = filename
        if(not pathlib.Path(template_file).exists()):
            print('Error: Template file {} not found.'.format(template_file))
            assert(False)
        return template_file

    def GetPythiaConfigFileContents(self,filename=None):
        filepath = self.GetPythiaConfigFile(filename)
        with open(filepath,'r') as f:
            contents = f.readlines()
        return ''.join(contents)

    # Get the Pythia configuration, for a single pT bin, as a dictionary.
    # This dictionary will only control a few settings (MPI, ISR/FSR etc.),
    # while the process settings from the 'proc' entry of config will be
    # passed separately.
    def GetPythiaConfig(self,pt_min, pt_max):
        pythia_config = {}
        # Tack on a few more things to the pythia configuration
        pythia_config['HadronLevel:all'] = Bool2String(self.config['hadronization'])
        pythia_config['PartonLevel:MPI'] = Bool2String(self.config['mpi'])
        pythia_config['PartonLevel:ISR'] = Bool2String(self.config['isr'])
        pythia_config['PartonLevel:FSR'] = Bool2String(self.config['fsr'])
        pythia_config['Random:setSeed'] = 'on'
        pythia_config['Random:seed'] = self.config['rng']

        # Add the phase space stuff.
        pythia_config['PhaseSpace:pTHatMin'] = pt_min
        pythia_config['PhaseSpace:pTHatMax'] = pt_max

        # Add some default (non-configurable) stuff.
        pythia_config['Print:quiet'] = 'on' # avoid printing reams of info
        pythia_config['Stat:showProcessLevel'] = 'off'
        pythia_config['Stat:showErrors'] = 'off'
    #     pythia_config['Next:numberShowProcess'] = '2'
    #     pythia_config['Next:numberShowEvent'] = '2'
        return pythia_config

    def GetPythiaRNGSeed(self):
        return self.config['rng']

    def GetTruthSelection(self):
        return self.config['truth_selection']

    def GetFinalStateSelection(self):
        return self.config['final_state_selection']

    def GetJetConfig(self):
        return_dict = {}
        for key in ['jet_radius','jet_min_pt','jet_max_eta','jet_n_par','jet_selection']: return_dict[key] = self.config[key]
        return return_dict

    def GetDelphesConfig(self):
        return self.config['delphes']

    def GetDelphesCard(self):
        return self.config['delphes_card']

    def GetNPars(self):
        return_dict = {}
        for key in ['jet_n_par','n_truth']: return_dict[key] = self.config[key]
        return return_dict

    def GetInvisiblesFlag(self):
        return self.config['invisibles']

    def GetEventSelection(self):
        return self.config['event_selection']

    def GetSignalFlag(self):
        return self.config['signal_flag']

    def GetSplitSeed(self):
        return self.config['split_seed']

    def GetRecordIndices(self):
        return self.config['record_final_state_indices']

    def GetPostProcessing(self):
        return self.config['post_processing']

    def GetEventFilter(self):
        return self.config['event_filter']

    def GetDelphesDirectory(self):
        return self.config['delphes_dir']

    def GetFastjetDirectory(self):
        return self.config['fastjet_dir']

    def GetPrintFastjet(self):
        return self.print_fastjet

    def GetPrintDelphes(self):
        return self.print_delphes

    def SetPrintFastjet(self,val):
        self.print_fastjet = val

    def SetPrintDelphes(self,val):
        self.print_delphes = val

    def GetUseVectorCalcs(self):
        return self.config['use_vectorcalcs']