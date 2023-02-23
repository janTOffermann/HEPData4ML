import sys,os
import pathlib
from config.config import config

def GetPythiaConfigFile():
    path_to_this = os.path.dirname(os.path.realpath(__file__))
    proc = config['proc']
    template_file = '{}/pythia_templates/{}.txt'.format(path_to_this,proc)
    if(not pathlib.Path(template_file).exists()):
        print('Error: Template file {} not found.'.format(template_file))
        assert(False)
    return template_file

def Bool2String(bool):
    if(bool): return 'on'
    return 'off'

# Get the Pythia configuration, for a single pT bin, as a dictionary.
# This dictionary will only control a few settings (MPI, ISR/FSR etc.),
# while the process settings from the 'proc' entry of config will be
# passed separately.
def GetPythiaConfig(pt_min, pt_max):

    pythia_config = {}

    # Tack on a few more things to the pythia configuration
    pythia_config['HadronLevel:all'] = Bool2String(config['hadronization'])
    pythia_config['PartonLevel:MPI'] = Bool2String(config['mpi'])
    pythia_config['PartonLevel:ISR'] = Bool2String(config['isr'])
    pythia_config['PartonLevel:FSR'] = Bool2String(config['fsr'])
    pythia_config['Random:setSeed'] = 'on'
    pythia_config['Random:seed'] = config['rng']

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

def GetTruthSelection():
    return config['truth_selection']

def GetFinalStateSelection():
    return config['final_state_selection']

def GetJetConfig():
    return_dict = {}
    for key in ['jet_radius','jet_min_pt','jet_max_eta','jet_n_par','jet_selection']: return_dict[key] = config[key]
    return return_dict

def GetDelphesConfig():
    return config['delphes']

def GetDelphesCard():
    return config['delphes_card']

def GetNPars():
    return_dict = {}
    for key in ['jet_n_par','n_truth']: return_dict[key] = config[key]
    return return_dict

def GetInvisiblesFlag():
    return config['invisibles']

def GetEventSelection():
    return config['event_selection']

def GetSignalFlag():
    return config['signal_flag']

def GetSplitSeed():
    return config['split_seed']

def GetRecordIndices():
    return config['record_final_state_indices']

def GetPostProcessing():
    return config['post_processing']