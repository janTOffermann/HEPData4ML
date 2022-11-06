import util.jet_selection as jetsel
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.event_selection as eventsel

config = {
    'proc' : 'Top_Wqq',
    'hadronization' : True,
    'mpi' : False,
    'isr' : False,
    'fsr' : False,
    'delphes' : True,
    'rng' : 3,
    'jet_radius': 0.8,
    'jet_min_pt': 15., #GeV
    'jet_max_eta': 2.,
    'jet_n_par': 200,
    'n_truth' : 155, # max number of truth particles to save per event
    # 'truth_selection' : truthsel.TruthSelection('t->Wb'),
    # 'truth_selection' : truthsel.AdvTruthSelection(parsel.SelectFinalStateDaughters(truthsel.TruthSelection('W')), 100),
    'truth_selection' : parsel.FirstSelector(1,22),
    # 'truth_selection' : None,
    'event_selection' : eventsel.TruthDistanceSelection(distance=2.4),
    # 'event_selection' : None,
    'final_state_selection': parsel.AlgoSelection(algos.SelectFinalState(),200),
    # 'final_state_selection' : parsel.SelectFinalStateDaughters(truthsel.TruthSelection('W')),
    'invisibles' : False, # if False, invisible particles in the final state will be discarded
    'jet_selection':jetsel.GetNearestJet(truth_code=6,max_dr=0.8),
    # 'jet_selection':jetsel.GetNearestJet(truth_code=24,max_dr=0.4),
    # 'jet_selection':jetsel.GetLeadingJet(),
    # 'jet_selection' : None,
    'signal_flag' : 1, # What to provide as the "signal_flag" for these events. (relevant if combining datasets). Must be >= 0.
    'split_seed' : 1 # seed to be used for the RNG when splitting dataset into train/test/validation samples
}

# Don't need to adjust the lines below.
try: config['truth_selection'].SetHadronization(config['hadronization'])
except: pass

# try: config['truth_selection'].SetN(config['n_truth'])
# except: pass
