# HEPData4ML

This package is meant to provide a relatively easy way to generate Monte Carlo datasets for high-energy physics (HEP) studies, such as (but not limited to) development of machine learning methods.

## Overview

This software package is designed for generating Monte Carlo simulated datasets of collider physics events.
It uses
- the [Pythia8](https://pythia.org) Monte Carlo generator to perform matrix element generation and showering/hadronization,
- the [DELPHES](https://cp3.irmp.ucl.ac.be/projects/delphes) library for fast detector simulation, and
- the [fastjet](http://fastjet.fr) library for jet clustering.
Future updates will allow for the supply of externally-produced HepMC3 files (as a replacement of the generation/showering step).

While the code on the main branch is currently structured to pick out a single jet per event, the `devel_full_event` branch features a much more general design, that saves entire events -- including possibly multiple jet collections that the user specifies.

Most options are specified in a Python configuration file (akin, for example, to the "steering" files used by the Athena software package in the ATLAS collaboration). By default, this resides in `config/config.py`, although alternative files can be used.

The final output of this package -- under normal running mode -- is an "n-tuple" in HDF5 format. Future updates may include support for a ROOT-based format, which may lend itself better to jagged arrays; the HDF5 files currently hold particle and jet properties in numpy arrays, which are truncated/zero-padded to user-specified sizes.

## Requirements/setup

The code is primarily written in Python3, along with some C++/ROOT code (which is leveraged in Python via PyROOT).
HEPData4ML makes use of a number of different software packages, including:

- [ROOT](https://root.cern.ch)
- [uproot v5](https://uproot.readthedocs.io/en/latest/index.html)
  - [awkward v2](https://awkward-array.org/doc/main/)
- [numpy](https://numpy.org)
- [h5py](https://www.h5py.org)
- [HepMC3](https://gitlab.cern.ch/hepmc/HepMC3)
- [Pythia8](https://pythia.org)
- [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes)
- [fastjet](http://fastjet.fr)

A YAML file detailing a `conda` environment with most of these dependencies is available in `setup/conda/`, including a binary distribution of ROOT.
For usage on CVMFS systems, the setup script in `setup/cvmfs/` is recommended instead, which will leverage `LCG_108`.
The Delphes, HepMC3 and fastjet dependencies will automatically be installed locally (i.e. within this repository) by the package during runtime, if they are not found (or their installation paths set to `None` in the Python configuration file); their locations should be provided in the aforementioned Python configuration file. A few notes:
- At the present moment, the HepMC3 path *should* be set to `None`, as we leverage features not yet available in an official HepMC3 release (but which are in its `master` branch on CERN GitLab).
- To make use of some custom Delphes

Note that the Pythia8 and fastjet installations that are used must have the Python interface built.

## Configuring and running the data generation

The data generation is run locally using the aptly-named `run.py` Python script. It takes a number of command line arguments:
```
usage: run.py [-h] -n NEVENTS [-steps STEPS] [-p [PTBINS ...]] [-o OUTFILE] [-O OUTDIR] [-v VERBOSE] [-f FORCE] [-c COMPRESS] [-rng RNG] [-npc NENTRIES_PER_CHUNK] [-pb PROGRESS_BAR] [-sp SPLIT] [-tf TRAIN_FRACTION] [-vf VAL_FRACTION]
              [-df DELETE_FULL] [-ds DELETE_STATS] [-co COMPRESSION_OPTS] [-pc PYTHIA_CONFIG] [-index_offset INDEX_OFFSET] [-config CONFIG] [-del_delphes DEL_DELPHES]

options:
  -h, --help            show this help message and exit
  -n, --nevents NEVENTS
                        Number of events per pt bin.
  -steps, --steps STEPS
                        Comma- or space-separated list of step. Options are [generation,simulation,reconstruction].
  -p, --ptbins [PTBINS ...]
                        Transverse momentum bin edges, for outgoing particles of the hard process. Can be a list of floats, or a string of comma- or space-separated floats. In GeV.
  -o, --outfile OUTFILE
                        Output HDF5 file name.
  -O, --outdir OUTDIR   Output directory.
  -v, --verbose VERBOSE
                        Verbosity.
  -f, --force FORCE     Whether or not to force generation -- if true, will possibly overwrite existing HepMC files in output directory.
  -c, --compress COMPRESS
                        Whether or not to compress HepMC files.
  -rng, --rng RNG       Pythia RNG seed. Will override the one provided in the config file.
  -npc, --nentries_per_chunk NENTRIES_PER_CHUNK
                        Number of entries to process per chunk, for jet clustering & conversion to HDF5.
  -pb, --progress_bar PROGRESS_BAR
                        Whether or not to print progress bar during event generation
  -sp, --split SPLIT    Whether or not to split HDF5 file into training/validation/testing files.
  -tf, --train_fraction TRAIN_FRACTION
                        Fraction of events to place in the training file.
  -vf, --val_fraction VAL_FRACTION
                        Fraction of events to place in the validation file.
  -df, --delete_full DELETE_FULL
                        Whether or not to delete the full HDF5 file after splitting into train/validation/testing files.
  -ds, --delete_stats DELETE_STATS
                        Whether or not to delete the full stats file. This file's info is merged into HDF5 dataset, but the file may be useful in some advanced use cases.
  -co, --compression_opts COMPRESSION_OPTS
                        Compression option for final HDF5 file (0-9). Higher value means more compression.
  -pc, --pythia_config PYTHIA_CONFIG
                        Path to Pythia configuration template (for setting the process).
  -index_offset, --index_offset INDEX_OFFSET
                        Offset for Event.Index.
  -config, --config CONFIG
                        Path to configuration Python file. Default will use config/config.py .
  -del_delphes, --del_delphes DEL_DELPHES
                        Whether or not to delete DELPHES/ROOT files.
```
A number of these are mostly likely not directly relevant to the user, and are mostly present for running using HTCondor (described further below). Note that many settings are provided via a Python configuration file (similar to how "steering file" work in ATLAS/Athena software, or CMSSW). For the command line arguments, here are descriptions of the most relevant ones:
- `--nevents`: How many events to produce, per pT bin (specified below).
- `--steps`: Here, one can list which particular steps in the event generation chain to run:
  - "generation" : This is done using Pythia8, and includes both matrix element generation and showering/hadronization. It produces HepMC3 files. This step can be skipped if there are existing HepMC3 files in the output directory, in which case those will be used in lieu of the generation output.
  - "pileup" : This step mixes pileup into the HepMC3 files, which is read from other (pre-existing) HepMC3 files. Note that these should specifically be HepMC3 files written in *ROOT* format, not ASCII (the generation step supports both output types). In principle, the pileup HepMC3 files can (and should) be produced first by a dedicated run of this software, using only the "generation" step. The output is HepMC3 files, with the pileup mixed into the event listings.
  - "simulation" : This step runs (fast) detector simulation, currently limited to Delphes. It outputs the "standard" Delphes output, consisitng ROOT files containing a TTree holding the Delphes n-tuple output. Note that technically some reconstruction is done in this step, despite the next step being called "reconstruction".
  - "reconstruction" : This step performs two tasks: it converts the output of previous steps (HepMC3 and optionally Delphes/ROOT) into HDF5, the final output format. It also runs any requested "post-processing" steps, which most importantly can include jet clustering (via the FastJet) library. The output is HDF5.
-  `--ptbins` : If running the "generation" step, this indicates the bins to use in "\hat{p}_{T}", the transverse momentum of the leading parton in event generation in Pythia8. A separate Pythia8 job is run per pT bin, in sequence.
- `--outfile` : The output file name. If the "split" option is used, this file name will be modified to contain suffixes "_train", "_test" and "_valid".
- `--outdir` : The output file directory.
- `--rng` : The random number generator (RNG) seed to use for generation with Pythia8. If the event generation is run twice with the same process and same RNG seed, it will produce the same events.
- `--split` : An integer (either `0` or `1`) indicating whether or not to split the output into training, testing and validation files. The fraction of events put into these can be managed via the `-train_fraction` and `--val_fraction` arguments.
- `--config` : A path to the Python configuration file to use, where most of the settings are actually provided. This defaults to `config/config.py`.
Note that certain arguments, such as `--pythia_config`, are redundant with the information provided in the Python configuration (described below). These arguments will override the settings in the Python configuration; their primary purpose is to facilitate overrides for HTCondor jobs (described further below).

### The Python configuration file
As noted above, most of the configuration is handled by a Python configuration file, by default that in `config/config.py`. Here is the annotated example file that ships with this repository:
```
import util.post_processing.jets as jets
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.pileup.pileup as pu

config = {
    'generation' : {
        'proc' : 'Top_Wqq', # Filename. Can also correspond to a card name in the util/pythia_templates subdirectory.
        'hadronization' : True, # Pythia8 hadronization flag
        'mpi' : True, # Pythia8 multi-parton interactions flag
        'isr' : True, # Pythia8 initial-state radiation flag
        'fsr' : True, # Pythia8 final-state radiation flag
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

    # NOTE: Object names (keys) should not have periods (".") in them. These are used internally to indicate objects' properties ("leaves" in ROOT-speak), and including these in names may break stuff down-the-line (such as in the visualization scripts).
    'reconstruction' : {
        'n_stable' : 200, # max number of stable truth-level particles to save per event (HDF5 doesn't support jagged arrays)
        'n_delphes': [200], # max number of Delphes objects to save per event -- list corresponding to entries in 'delphes_output'. If single value, will be broadcast to appropriate shape.
        'fastjet_dir' : None, # Directory containing the Fastjet installation. If None, will build in a local directory "external/fastjet". Note that the Fastjet installation you use must have the Python bindings set up.
        'n_truth' : 1 + 60, # Maximum number of truth particles to save per event. (HDF5 doesn't support jagged arrays)
        'event_filter' : None, # Deprecated.
        'event_filter_flag': None, # Deprecated.
        'particle_selection' : { # Here, you can specify collections of truth-level particles to save, using various (provided) particle selection algorithms. These typically search for particles matching some PdgID and/or generator status.
            'TruthParticlesTopAndChildren':
            parsel.MultiSelection(
                [
                    parsel.FirstSelector(22, 6), # top quark
                    parsel.FirstSelector(23, 5), # bottom quark
                    parsel.FirstSelector(22,24), # W boson
                    parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(22,24)),n=120) # up to 120 stable daughters of W
                ]
            ),
            'TruthParticlesAntiTopAndChildren':
            parsel.MultiSelection(
                [
                    parsel.FirstSelector(22, -6), # top anti-quark
                    parsel.FirstSelector(23, -5), # bottom anti-quark
                    parsel.FirstSelector(22,-24), # W boson
                    parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(22,-24)),n=120) # up to 120 stable daughters of W
                ]
            )
        },
        'signal_flag' : 1, # What to provide as the "SignalFlag" for these events. Relevant if combining multiple samples, so as to bookkeep what is what.
        'split_seed' : 1, # RNG seed to be used for splitting the dataset into train/test/validation samples.
        'post_processing': [ # What post-processing algorithms to run -- this includes jet clustering! You can queue up multiple separate post-processors.

            # Cluster large-radius jets, ghost associated to the top and antitop. Enforce some pt and eta cuts.
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08RecoJetsAssociatedTop').PtFilter(25.).EtaFilter(4.).GhostAssociation('TruthParticlesTopAndChildren',    0,mode='filter'),
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08RecoJetsAssociatedAntiTop').PtFilter(25.).EtaFilter(4.).GhostAssociation('TruthParticlesAntiTopAndChildren',0,mode='filter'),

            # Cluster small-radius jets, near the ghost-associated jets above.
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04RecoJetsAssociatedTop').Containment('AntiKt08RecoJetsAssociatedTop',0,0.4,mode='filter'),
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.4,jet_name='AntiKt04RecoJetsAssociatedAntiTop').Containment('AntiKt08RecoJetsAssociatedAntiTop',0,0.4,mode='filter')
        ]
    }
}
```
As you can see, the general format is a Python dictionary, with sub-dictionaries associated with each step in the event generation chain.

### Running with HTCondor
In practice, you will probably want to run a large number of event generation jobs in parallel, as this is inherently a very parallelizable task. This can be accomplished using the `prep_condor.py` and `run_condor.py` scripts; the former prepares a condor submission script with the accompanying arguments, and the latter simply submits that set of jobs (and generates the necessary job subdirectories for log files). Here are the arguments that `prep_condor.py` takes:
```
usage: prep_condor.py [-h] -n NEVENTS [-steps STEPS] [-p [PTBINS ...]] [-O OUTDIR] [-R RUNDIR] [-rng RNG] [-sp SPLIT] [-pc PYTHIA_CONFIG [PYTHIA_CONFIG ...]] [-N NJOBS] [-sq SHORTQUEUE] [-ncpu N_CPU] [-nblas N_OPENBLAS]
                      [-event_idx_offset EVENT_IDX_OFFSET] [-batch_name BATCH_NAME] [-mode MODE] [-branch BRANCH] [-config CONFIG] [-template TEMPLATE] [-requirements REQUIREMENTS] [-blacklist BLACKLIST]

options:
  -h, --help            show this help message and exit
  -n, --nevents NEVENTS
                        Number of events per pt bin.
  -steps, --steps STEPS
                        Comma- or space-separated list of step. Options are [generation,simulation,reconstruction].
  -p, --ptbins [PTBINS ...]
                        Transverse momentum bin edges, for outgoing particles of the hard process. Can be a list of floats, or a string of comma- or space-separated floats. In GeV.
  -O, --outdir OUTDIR   Output directory for the jobs
  -R, --rundir RUNDIR   Run directory -- where the condor jobs will go.
  -rng, --rng RNG       Pythia RNG seed offset -- seed will be offset + job number.
  -sp, --split SPLIT    Whether or not to split HDF5 file into training/validation/testing files.
  -pc, --pythia_config PYTHIA_CONFIG [PYTHIA_CONFIG ...]
                        Path to Pythia configuration template (for setting the process).
  -N, --Njobs NJOBS     Number of jobs per config.
  -sq, --shortqueue SHORTQUEUE
                        Whether or not to use the condor short queue (for UChicago Analysis Facility).
  -ncpu, --n_cpu N_CPU  Number of (logical) CPU cores per job.
  -nblas, --n_openblas N_OPENBLAS
                        Sets the $OPENBLAS_NUM_THREADS variable (used for numpy multithreading). Advanced usage.
  -event_idx_offset, --event_idx_offset EVENT_IDX_OFFSET
                        Initial offset for event_idx. Advanced usage.
  -batch_name, --batch_name BATCH_NAME
                        Name for the condor batch. If left empty, unused.
  -mode, --mode MODE    If <=0, ships code directly from this repo. If 1, runs 'git clone' within each job. If 2, runs 'git clone' locally and ships result to jobs. If 3, points jobs to this repo (doesn't ship it).
  -branch, --branch BRANCH
                        If using git for jobs, which branch to clone. Defaults to current.
  -config, --config CONFIG
                        Path to Python config file, for run.py .
  -template, --template TEMPLATE
                        Template condor submission file. Templates are in util/condor_templates.
  -requirements, --requirements REQUIREMENTS
                        Requirements string for condor workers.
  -blacklist, --blacklist BLACKLIST
                        Blacklist file for condor workers.
```
Some of these arguments directly match those given to `run.py`. There are a few that are however specific to this script, and have to do with how the condor jobs are set up:
- `--Njobs`: The number of jobs to run, per Pythia8 configuration.
- `--rundir`: The directory from which the jobs will be run (more on this below, in the description of the `--mode` argument). This is where the condor submission file will be created, as well as where the job subdirectories that contain condor logs and stdout/stderr will be produced when `run_condor.py` is invoked.
- `--pythia_config`: A path to the Pythia8 configuration file; this overrides whatever configuration is set in the Python configuration file, however this can be a *list* of Pythia8 configuration files in which case `Njobs` jobs will be spun up per Pythia8 configuration. This can be a handy way to simultaneously set up generation of multiple processes. Note however that this may make it difficult to later distinguish between these processes, whereas spinning up separate batches of jobs per process will allow for labeling them via the `signal_flag` setting in the `reconstruction` section of the Python configuration file.
- `--shortqueue` : This setting is currently unique to the HTCondor system on the University of Chicago's Analysis Facility, where it allows for accessing a dedicated HTCondor queue for jobs limited to 3 hours. (This code can possibly be extended to similarly access unique features of other HTCondor systems).
- `--n_cpu`: The number of CPU cores to use per job. Note that this code is not (yet) explicitly parallelized, so there is little benefit to requesting multiple cores (aside from that possibly also providing more memory).
- `batch_name`: The name for this batch of condor jobs.
- `mode`: Which "run mode" to use; this is a consequence of HTCondor queues operating in slightly different ways on different clusters. Here are the available options:
  - `1`: The condor jobs will each clone the `HEPData4ML` repository. This requires them to have the ability to clone repositories from GitHub, and any uncommitted local changes will not be present (the clone will use whatever branch of the code you are currently sitting on).
  - `2`: A git clone of the repository will be performed locally, and then this will be shipped to the condor jobs  as a `tar.gz` archive using HTCondor's file transfer protocol.
  - `3`: The code will run "locally", i.e. the condor workers will run the code directly from this repository. This requires the workers to have read access to wherever this repository is located. This mode is useful for example on Brown University's BRUX cluster,  where the HTCondor jobs are not "shipped off" to scratch directories by default.
  - For any other value, the local repository will be packaged into a `tar.gz` archive, and shipped to the condor jobs usingHTCondor's file transfer protocol. This is the "typical" run mode.
- `requirements`: A string of requirements for the condor worker nodes.
- `blacklist`: Specific condor worker machines to avoid running on; this can be a handy feature for systems where certain worker nodes don't function well, e.g. some nodes not having consistent access to CVMFS.

**Note** that at the present, the condor workers must have write access to the output directory. This is something that may be changed in the future.

### Parsing and handling output
The final output format of the data generation is an HDF5 file, that contains an n-tuple representing the data. It comprises of a set of n-dimensional `numpy` arrays, which represent properties of the event; the four-momenta of truth-level particles, the four-momenta (and position, charge etc.) of various reconstructed objects from Delphes as well as jet clustering.

HDF5 has the advantage that it is quite widely supported (and relatively popular within the ML community). However unlike ROOT it does not support jagged arrays, which are a relatively natural representation for a lot of HEP data as (for example) the number of particles per event, or number of constituents per jet, is a variable quantity. Thus the n-tuples are constructed using fixed-length arrays, where the array length is (hopefully!) chosen to be sufficiently large as to avoid any edge effects, and zero-padding is employed as needed; each object in the n-tuple also has a column indicating its "multiplicity", i.e. how many of them there are per event.

In addition to the kinematics and properties of objects in each event, there are also fields corresponding with *metadata*, such as the configuration used to generate the events, with these fields indexing each event with respect to a list of metadata values stored in the HDF5 file's "attributes". When files are concatenated using the Python script `util/tools/concat.py`, merging of the metadata is handled automatically. As it includes the configuration, command line arguments and git hash (among other pieces of information), this metadata in principle allows one to completely reproduce a dataset (aside from the Delphes RNG -- which may be handled in future updates!).

One can check the particular contents of a dataset by using the Python script `util/tools/check_file.py`, which will print all the available fields, their dimensions, as well as a single event.

