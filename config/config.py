import numpy as np
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.event_filter as eventfilter
# import util.post_processing as postproc
import util.post_processing.jets as jets
# import util.post_processing.tracing as tracing
# import util.post_processing.truth_sum as truth_sum
# import util.post_processing.rotate as rotate
# import util.post_processing.containment as containment
# import util.post_processing.jh_tagging as jh_tagging

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
    'particle_selection' : {
        'TruthParticles':
        parsel.MultiSelection(
            [
                parsel.FirstSelector(23,5), # bottom quark
                parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(23,5)),n=60) # up to 60 stable daughters of b quark
            ]
        ),
    },
    'signal_flag' : 1, # What to provide as the "signal_flag" for these events. (relevant if combining datasets). Must be >= 0.
    'split_seed' : 1, # seed to be used for the RNG when splitting dataset into train/test/validation samples
    # 'use_vectorcalcs' : False, # whether or not to use the VectorCalcs C++/ROOT library (which is part of this repo). May speed up some computations, but can lead to issues if there's a problem with the ROOT build (e.g. CVMFS/LCG_103 seems to cause issues)
    'post_processing': [
        # cluster jets
        jets.JetFinder('StableTruthParticles',jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04GenJets'),
        jets.JetFinder('StableTruthParticles',jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04GenJets_bGhostAssociated').GhostAssociation('TruthParticles',0),
        jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04EFlowJets'),
        jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08EFlowJets')
    ]
}

# Don't adjust the lines below.
try: config['truth_selection'].SetHadronization(config['hadronization'])
except: pass
