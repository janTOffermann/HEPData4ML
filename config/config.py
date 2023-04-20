import util.jet_selection as jetsel
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.event_selection as eventsel
import util.event_filter as eventfilter
import util.post_processing.tracing as tracing
import config.selectors as s

selections = s.selections

config = {
    'proc' : 'Top_Wqq',
    'hadronization' : True,
    'mpi' : False,
    'isr' : False,
    'fsr' : False,
    'rng' : 1, # Pythia RNG seed
    'delphes' : False, # Whether or not to use Delphes
    'delphes_card' : None, # path to the Delphes card to use. If None, will use the ATLAS Delphes card that ships with Delphes
    # 'delphes_dir' : None ,# directory containing the Delphes installation. If None, will be build in a local directory "delphes".
    'delphes_dir' : "/cvmfs/sft.cern.ch/lcg/releases/delphes/3.5.1pre05-775ca/x86_64-centos7-gcc11-opt",
    'fastjet_dir' : "/home/jaofferm/HEPData4ML/fastjet",
    #'fastjet_dir' : None,
    'jet_radius': 0.8,
    'jet_min_pt': 15., #GeV # 15
    'jet_max_eta': 2., # absolute value eta cut # 2.
    'jet_n_par': 200, # max number of jet constituents to save per jet
    'n_truth' : 3 + 2 + 120, # max number of truth particles to save per jet
    # 'event_filter' : None,
    'event_filter' : eventfilter.MultiFilter(
        [
            eventfilter.PtMatchedJetFilter(0.8,selections['t->Wb'],pt_window_frac=0.005,pt_min_jet=15.,eta_max_jet=2.),
            eventfilter.ContainedJetFilter(0.8,selections['t->Wb'],selections['bqq'],matching_radius=0.6,pt_min_jet=15.,eta_max_jet=2.),
            eventfilter.ContainedJetFilter(0.8,selections['t->Wb'],selections['t daughters'],pt_threshold_frac=0.01,pt_min_jet=15.,eta_max_jet=2.)
        ]
    ),
    'truth_selection' : selections['t->Wb w/ qq and W daughters'],
    'final_state_selection': parsel.AlgoSelection(algos.SelectFinalState(),200),
    'event_selection' : eventsel.TruthDistanceSelection(distance=2.4, n_truth=3), # filters an event to remove some final-state particles, primarily for lowering HepMC file-size. Ideally does not affect the final output.
    'invisibles' : False, # if False, invisible particles in the final-state selection (neutrinos) will be discarded. If true, they may be input to jet clustering!
    'jet_selection':jetsel.GetNearestJet(truth_code=6,max_dr=0.8),
    # 'jet_selection':jetsel.GetNearestJetWithContainment(6,0.8,5,-1),
    'signal_flag' : 1, # What to provide as the "signal_flag" for these events. (relevant if combining datasets). Must be >= 0.
    'post_processing': [
        None,
        # tracing.Tracer(verbose=False)
    ],
    'record_final_state_indices' : True, # Whether or not to record jet constituents' indices w.r.t. the order they were passed to jet clustering (order of particles in HepMC file, or order of Delphes objects if using Delphes).
    'split_seed' : 1, # seed to be used for the RNG when splitting dataset into train/test/validation samples
    'use_vectorcalcs' : False # whether or not to use the VectorCalcs C++/ROOT library (which is part of this repo). May speed up some computations, but can lead to issues if there's a problem with the ROOT build (e.g. CVMFS/LCG_103 seems to cause issues)
}

# Don't need to adjust the lines below.
try: config['truth_selection'].SetHadronization(config['hadronization'])
except: pass
