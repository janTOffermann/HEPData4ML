import util.post_processing.jets as jets
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.pileup.pileup as pu

config = {
    'generation' : {
        'proc' : 'SoftQCD', # Filename. Can also correspond to a card name in the util/pythia_templates subdirectory.
        'hadronization' : True, # Pythia8 hadronization flag
        'mpi' : True, # Pythia8 multi-parton interactions flag
        'isr' : True, # Pythia8 initial-state radiation flag
        'fsr' : True, # Pythia8 final-state radiation flag
        'rng' : 1, # Pythia8 RNG seed
        'verbose' : False,
        'hepmc_dir': None, # Directory containing the HepMC3 installation. If None, will build in a local directory "external/hepmc". Note that we currently use a custom fork of HepMC3, so you probably want to leave this as None.
        'hepmc_format': 'root' # Options are 'root' and 'ascii'. The ROOT option provides superior filesize and random-access capability (useful if making samples for pileup overlay), at the cost of being less human-readable.
    },

    # When using this config file, we'll be sure to just do the "generation" step, so we can safely omit
    # the sections for "pileup", "simulation" and "reconstruction".
}
