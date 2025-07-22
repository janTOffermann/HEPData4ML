import util.post_processing.jets as jets
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.pileup.pileup as pu

config = {
    'generation' : {
        'proc' : 'Top_Wqq', # Filename. Can also correspond to a card name in the util/pythia_templates subdirectory.
        'hadronization' : True, # Pythia8 hadronization flag
        'mpi' : False, # Pythia8 multi-parton interactions flag
        'isr' : False, # Pythia8 initial-state radiation flag
        'fsr' : False, # Pythia8 final-state radiation flag
        'rng' : 1, # Pythia8 RNG seed
        'verbose' : False,
        'hepmc_dir': None, # Directory containing the HepMC3 installation. If None, will build in a local directory "external/hepmc". Note that we currently use a custom fork of HepMC3, so you probably want to leave this as None.
        'hepmc_format': 'root' # Options are 'root' and 'ascii'. The ROOT option provides superior filesize and random-access capability (useful if making samples for pileup overlay), at the cost of being less human-readable.
    },

    'pileup' : {
        'handler':None,
        # 'handler': pu.PileupOverlay("/Users/jan/tmp/pileup/part0/events_0.root",rng_seed=1) # For example, you can overlay pileup events from some pre-existing HepMC3 files (ideally in ROOT format!), which you can generate with this package too.
    },

    'simulation' : {
        'type' : 'delphes', # what simulation (if any) to use. Currently supported options are [None, 'delphes']
        'delphes_card' : "util/delphes/cards/delphes_card_CMS_custom.tcl", # path to the Delphes card to use. If None, will use the ATLAS Delphes card that ships with Delphes
        'delphes_dir' : None, # Directory containing the Delphes installation. If None, will be build in a local directory "external/delphes". If using Delphes from CVMFS, this should match your release/views setup, otherwise it might not work! E.g. "/cvmfs/sft.cern.ch/lcg/releases/LCG_105/delphes/3.5.1pre09/x86_64-el9-gcc13-opt". Note that our custom CMS card requires a custom fork of Delphes (which will be installed if None).
        'delphes_output' : ['EFlowPhoton','EFlowNeutralHadron','EFlowTrack','Electron','Muon','Photon','GenMissingET','MissingET','GenVertex','Vertex'] # Which output objects from Delphes to propagate to the final HDF5 file -- this is also what will be available to the post-processors; other information will be dropped. Some details for the vertex-type objects may still need some ironing out.
    },

    'reconstruction' : {
        'n_stable' : 200, # max number of stable truth-level particles to save per event (HDF5 doesn't support jagged arrays)
        'n_delphes': [200], # max number of Delphes objects to save per event -- list corresponding to entries in 'delphes_output'. If single value, will be broadcast to appropriate shape.
        'fastjet_dir' : None, # Directory containing the Fastjet installation. If None, will build in a local directory "external/fastjet". Note that the Fastjet installation you use must have the Python bindings set up.
        'n_truth' : 1 + 60, # Maximum number of truth particles to save per event. (HDF5 doesn't support jagged arrays)
        'event_filter' : None, # Deprecated.
        'event_filter_flag': None, # Deprecated.
        'particle_selection' : { # Here, you can specify collections of truth-level particles to save, using various (provided) particle selection algorithms. These typically search for particles matching some PdgID and/or generator status.
            'TruthParticles':
            parsel.MultiSelection(
                [
                    parsel.FirstSelector(22, 6), # top quark
                    parsel.FirstSelector(23, 5), # bottom quark
                    parsel.FirstSelector(22,24), # W boson
                    parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(22,24)),n=120) # up to 120 stable daughters of W
                ]
            ),
        },
        'signal_flag' : 1, # What to provide as the "SignalFlag" for these events. Relevant if combining multiple samples, so as to bookkeep what is what.
        'split_seed' : 1, # RNG seed to be used for splitting the dataset into train/test/validation samples.
        'post_processing': [ # What post-processing algorithms to run -- this includes jet clustering! You can queue up multiple separate post-processors.
            # cluster jets
            # jets.JetFinder('StableTruthParticles',jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04GenJets').PtFilter(15.), # some truth jets, with at least 15 GeV in pT
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04RecoJets').PtFilter(15.).GhostAssociation('TruthParticles',0,mode='tag'), # some reco-level jets -- incl. a tag for whether or not they're ghost-associated with the bottom quark selected above
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08RecoJets').PtFilter(15.), # some reco-level jets -- larger radius

            # jets.JetFinder('StableTruthParticles',jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08GenJets').PtFilter(15.).EtaFilter(2.).JohnsHopkinsTagger(mode='tag')
            # jets.JetFinder('StableTruthParticles',jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04GenJets_bGhostAssociated').GhostAssociation('TruthParticles',0,mode='tag'),
            # jets.JetFinder('StableTruthParticles',jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08GenJets').JohnsHopkinsTagger(mode='tag')
            # jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04EFlowJets'),
            # jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08EFlowJets')
        ]
    }
}
