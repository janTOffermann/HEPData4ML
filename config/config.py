import numpy as np
import util.jet_selection as jetsel
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.event_selection as eventsel
import util.event_filter as eventfilter
import util.post_processing.tracing as tracing
import util.post_processing.truth_sum as truth_sum
import util.post_processing.rotate as rotate
import util.post_processing.containment as containment
import util.post_processing.jh_tagging as jh_tagging
import config.selectors as s

selections = s.selections

config = {
    'proc' : 'bbar', # Name must correspond to a card name (without extension) in util/pythia_templates
    'hadronization' : True, # Pythia8 hadronization flag
    'mpi' : False, # Pythia8 multi-parton interactions flag
    'isr' : False, # Pythia8 initial-state radiation flag
    'fsr' : False, # Pythia8 final-state radiation flag
    'rng' : 1, # Pythia8 RNG seed
    'delphes' : True, # Whether or not to use Delphes
    'delphes_card' : "util/delphes/cards/delphes_card_CMS.tcl", # path to the Delphes card to use. If None, will use the ATLAS Delphes card that ships with Delphes
    'delphes_dir' : '/home/jofferma/projects/HEPData4ML/delphes',# '/cvmfs/sft.cern.ch/lcg/releases/LCG_105/delphes/3.5.1pre09/x86_64-el9-gcc13-opt' ,# directory containing the Delphes installation. If None, will be build in a local directory "delphes". For CVMFS, this should match your release/views, otherwise might not work!
    'delphes_output' : ['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],
    'n_stable' : 200, # max number of stable truth-level particles to save to HDF5 file
    'n_delphes': [200], # max number of Delphes objects to save to HDF5 file -- list corresponding to entries in 'delphes_output'
    'fastjet_dir' : '/home/jofferma/projects/HEPData4ML/fastjet',
    'n_truth' : 1 + 60, # max number of truth particles to save per jet
    'event_filter' : None, # if not using any event filter
    'event_filter_flag': None, # an event_filter algorithm, but instead of filtering out events it will simply add a boolean flag to the dataset
    'truth_selection' : parsel.MultiSelection(
        [
            parsel.FirstSelector(23,5), # bottom quark
            # parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(23,5)),n=60) # up to 60 stable daughters of b quark
        ]
    ),
    'final_state_selection': parsel.AlgoSelection(algos.SelectFinalState(),200),
    'invisibles' : False, # if False, invisible particles in the final-state selection (neutrinos) will be discarded. If true, they may be input to jet clustering (which may be an issue at truth-level, though not DELPHES)!
    'jet_selection':jetsel.GetNearestJet(truth_code=5,max_dr=0.4), # pick jet nearest truth-level b, and within dR=0.4
    'signal_flag' : 1, # What to provide as the "signal_flag" for these events. (relevant if combining datasets). Must be >= 0.
    'record_final_state_indices' : True, # Whether or not to record jet constituents' indices w.r.t. the order they were passed to jet clustering (order of particles in HepMC file, or order of Delphes objects if using Delphes).
    'split_seed' : 1, # seed to be used for the RNG when splitting dataset into train/test/validation samples
    'use_vectorcalcs' : False, # whether or not to use the VectorCalcs C++/ROOT library (which is part of this repo). May speed up some computations, but can lead to issues if there's a problem with the ROOT build (e.g. CVMFS/LCG_103 seems to cause issues)
    'post_processing': None
}

# Don't adjust the lines below.
try: config['truth_selection'].SetHadronization(config['hadronization'])
except: pass
