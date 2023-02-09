import util.jet_selection as jetsel
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.event_selection as eventsel
import config.selectors as s

selections = s.selections

config = {
    'proc' : 'Top_Wqq',
    'hadronization' : True,
    'mpi' : False,
    'isr' : False,
    'fsr' : False,
    'delphes' : True,
    'rng' : 1, # Pythia RNG seed
    'jet_radius': 0.8,
    'jet_min_pt': 15., #GeV
    'jet_max_eta': 2., # absolute value eta cut
    'jet_n_par': 200, # max number of jet constituents to save per jet
    'n_truth' : 3 + 2 + 200, # max number of truth particles to save per jet
    'truth_selection' : selections['t->Wb w/ qq and daughters'],
    'final_state_selection': parsel.AlgoSelection(algos.SelectFinalState(),200),
    'event_selection' : eventsel.TruthDistanceSelection(distance=2.4, n_truth=3), # filters an event to remove some final-state particles, primarily for lowering HepMC file-size. Ideally does not affect the final output.
    'invisibles' : False, # if False, invisible particles in the final state will be discarded
    'jet_selection':jetsel.GetNearestJet(truth_code=6,max_dr=0.8),
    'signal_flag' : 1, # What to provide as the "signal_flag" for these events. (relevant if combining datasets). Must be >= 0.
    'record_final_state_indices' : True, # Whether or not to record jet constituents' indices w.r.t. the order they were passed to jet clustering (order of particles in HepMC file, or order of Delphes objects if using Delphes).
    'split_seed' : 1 # seed to be used for the RNG when splitting dataset into train/test/validation samples
}

# Don't need to adjust the lines below.
try: config['truth_selection'].SetHadronization(config['hadronization'])
except: pass
