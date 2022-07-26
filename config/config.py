import util.jet_selection as jetsel
import util.truth_selection as truthsel
import util.particle_selection as parsel

config = {
    'proc': 'Top_Wqq',
    'hadronization':False,
    'mpi':False,
    'isr':False,
    'fsr':False,
    'delphes':False,
    'rng':1,
    'jet_radius':0.8,
    'jet_min_pt':15., #GeV
    'jet_max_eta':2.,
    'jet_n_par': 3, # 200
    'n_truth':3,
    'truth_selection': truthsel.TruthSelection('t->Wb_nohad'),
    # 'final_state_selection': parsel.SelectFinalState,
    'final_state_selection': parsel.SelectSimplestHadronic(truthsel.TruthSelection('Wb_nohad')),
    # 'jet_selection':jetsel.GetTopJet
    'jet_selection':None,
    'alpha':3. # for filtering HepMC files -- particles this many jet radii away from all selected truth particles will be instantly discarded. Should be >> 1.
}

# Don't need to adjust the lines below.
config['truth_selection'].SetHadronization(config['hadronization'])