import util.jet_selection as jetsel
import util.truth_selection as truthsel
import util.particle_selection as parsel
import util.event_selection as eventsel

config = {
    'proc' : 'Top_Wqq',
    'hadronization' : True,
    'mpi' : False,
    'isr' : False,
    'fsr' : False,
    'delphes' : False,
    'rng' : 1,
    'jet_radius': 0.8,
    'jet_min_pt': 15., #GeV
    'jet_max_eta': 2.,
    'jet_n_par': 200,
    'n_truth' : 3,
    'truth_selection' : truthsel.TruthSelection('t->Wb'),
    'event_selection' : eventsel.TruthDistanceSelection(distance=2.4),
    'final_state_selection': parsel.SelectFinalState,
    'invisibles' : False, # if False, invisible particles in the final state will be discarded
    # 'final_state_selection' : parsel.SelectSimplestHadronic(truthsel.TruthSelection('Wb_nohad')),
    'jet_selection':jetsel.GetNearestJet(truth_code=6,max_dr=0.8),
    # 'jet_selection' : None,
    # 'event_selection' : None
}

# Don't need to adjust the lines below.
config['truth_selection'].SetHadronization(config['hadronization'])