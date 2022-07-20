import sys,os
import pathlib
#import numpythia as npyth # for defining selections
import util.jet_selection as jetsel
import util.truth_selection as truthsel

config = {
    'proc': 'Top_Wqq',
    'hadronization':True,
    'mpi':False,
    'isr':False,
    'fsr':False,
    'delphes':False,
    'rng':1,
    'jet_radius':0.8,
    'jet_min_pt':15., #GeV
    'jet_max_eta':2.,
    'jet_n_par': 200,
    'n_truth':3,
    'truth_selection': truthsel.selections['t->Wb'],
    'jet_selection':jetsel.GetTopJet
}

def GetPythiaConfigFile():
    global config
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
    global config

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

def GetJetConfig():
    return_dict = {}
    for key in ['jet_radius','jet_min_pt','jet_max_eta','jet_n_par','jet_selection']: return_dict[key] = config[key]
    return return_dict

def GetDelphesConfig():
    return config['delphes']

def GetNPars():
    return_dict = {}
    for key in ['jet_n_par','n_truth']: return_dict[key] = config[key]
    return return_dict