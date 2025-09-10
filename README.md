# HEPData4ML

This package is meant to provide a way to generate Monte Carlo datasets for high-energy physics (HEP) studies, such as (but not limited to) development of machine learning methods.

Plenty of solutions to this problem already exist; this package is meant to provide a relatively easy-to-use one, and produces output files replete with plenty of metadata to help document how they were created (something all too often missing from public benchmark datasets!).

## Overview

This software package is designed for generating Monte Carlo simulated datasets of collider physics events.
It uses
- the [Pythia8](https://pythia.org) Monte Carlo generator to perform matrix element generation and showering/hadronization,
- the [DELPHES](https://cp3.irmp.ucl.ac.be/projects/delphes) library for fast detector simulation, and
- the [fastjet](http://fastjet.fr) library for jet clustering.
The package is designed for flexibility and steps can be excluded from the workflow; e.g. one can start with pre-existing HepMC3 files, and just run fast detector simulation and object reconstruction.

Much of the configuration is specified in a Python configuration file (herein referred to as the "config file") akin, for example, to the steering files used by the Athena software package in the ATLAS collaboration. By default, this resides in `config/config.py`, although alternative files can be used. Being a Python file, the config file can hold complex objects -- and thus is a nice way to define what can be pretty intricate configurations for how what types of particles you want to save to output, how to reconstruct jets, et cetera.

The final output of this package -- under normal running mode -- is an "n-tuple" in HDF5 format. Future updates may include support for a ROOT-based format, which may lend itself better to jagged arrays; the HDF5 files currently hold particle and jet properties in numpy arrays, which are truncated/zero-padded to user-specified sizes.

## Requirements/setup

This package makes use of Git submodules; after cloning, please run
```
git submodule update --init --recursive
```
to make sure that these are fetched as well.


### Software requirements overview
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

The software environment can be set up using the `conda` package manager, or by using pre-packaged environments available via the CERN Virtual Machine File System (CVMFS).

### Conda setup

A YAML file detailing a `conda` environment with most of these dependencies is available in `setup/conda/`, including a binary distribution of ROOT. This environment can be set up by running
```
conda create -f setup/conda/hepdata4ml.yml # check this!
```
### CVMFS usage
For usage on CVMFS systems, the setup script in `setup/cvmfs/` is recommended instead, which will leverage `LCG_108`.
```
source setup/cvmfs/setup.sh
```

### Notes on dependencies: automatically installed packages and git submodules
The Delphes, HepMC3 and fastjet dependencies will automatically be installed in `external/` (i.e. within this repository) by the package during runtime, if they are not found (or their installation paths set to `None` in the config file), with the former two being `git submodules` of this package.
Should you wish to intead use existing installations of these packages, their locations should be provided in the aforementioned config file. A few notes:
- The code leverages features in HepMC3 and Delphes that may not presently be available in official releases (or in the case of Delphes, in the official repository: the submodule corresponds to my own fork of the project).
    - The "non-standard" HepMC3 features are presently part of the master branch of the [official HepMC3 repository](https://gitlab.cern.ch/hepmc/HepMC3/-/tree/master?ref_type=heads), so it is foreseeable that an official release may eventually be used instead.
    - By contrast, the Delphes-based code makes use of some special features that are only available in my fork of the project; this is true both of some of the provided Delphes detector cards, as well as the `DelphesHepMC3ROOT` executable in `external/delphes_custom` (as well as the `DelphesWrapper` Python interface).

Note that the Pythia8 and fastjet installations that are used *must* have the Python interface built -- running the script will check the fastjet Python interface, and if it's not present it will install in `external/fastjet`.

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
Note that a number of these command-line arguments are likely not directly relevant to the user; they are there for internal use when producing HTCondor jobs to run in parallel (described further below). In addition, some of these arguments are redundant with those in the config file: if provided, they will override the corresponding config file contents. Here are descriptions of the most relevant ones:
- `--nevents`: How many events to produce, per pT bin (specified below).
- `--steps`: Here, one can list which particular steps in the event generation chain to run:
  - "generation" : This is done using Pythia8, and includes both matrix element generation and showering/hadronization. It produces HepMC3 files. This step can be skipped if there are existing HepMC3 files in the output directory, in which case those will be used in lieu of the generation output.
  - "pileup" : This step mixes pileup into the HepMC3 files, which is read from other (pre-existing) HepMC3 files. Note that these should specifically be HepMC3 files written in *ROOT* format, not ASCII (the generation step supports both output types). In principle, the pileup HepMC3 files can (and should) be produced first by a dedicated run of this software, using only the "generation" step. The output is HepMC3 files, with the pileup mixed into the event listings.
  - "simulation" : This step runs (fast) detector simulation, currently limited to Delphes. It outputs the "standard" Delphes output, consisting of ROOT files containing a TTree holding the Delphes n-tuple output. Note that technically some reconstruction is done in this step, despite the next step being called "reconstruction".
  - "reconstruction" : This step performs two tasks: it converts the output of previous steps (HepMC3 and optionally Delphes/ROOT) into HDF5, the final output format. It also runs any requested "post-processing" steps, which most importantly can include jet clustering (via the FastJet) library. The output is HDF5.
-  `--ptbins` : If running the "generation" step, this indicates the bins to use in "\hat{p}_{T}", the transverse momentum of the leading parton in event generation in Pythia8. A separate Pythia8 job is run per pT bin, in sequence.
- `--outfile` : The output file name. If the "split" option is used, this file name will be modified to contain suffixes "_train", "_test" and "_valid".
- `--outdir` : The output file directory.
- `--rng` : The random number generator (RNG) seed to use for Pythia8, pileup mixing and Delphes (if relevant). Note that in general, if the event generation is run twice with the same settings and the same RNG seed, it will produce the same events: modifying the RNG seed is an easy way to guarantee that one will produce different events, and thus is a good way of "extending" a dataset.
- `--split` : If provided, a flag indicating whether or not to split the output into training, testing and validation files. The fraction of events put into these can be managed via the `-train_fraction` and `--val_fraction` arguments, which are floats on the interval `[0.,1.)`. Events that don't go into the training or validation splits will be placed in the testing split.
- `--config` : A path to the config file to use, where most of the settings are actually provided. This defaults to `config/config.py`.

### The config file
As noted above, most of the configuration is handled by a config file, by default that in `config/config.py`. Here is the annotated example file that ships with this repository:
```
import util.reconstruction.post_processing.jets as jets
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.pileup.pileup as pu

# The configuration is provided as a dictionary of dictionaries. The "sub-dictionaries" correspond with the different available steps:
# ['generation','pileup','simulation','reconstruction'].
config = {
    'generation' : { # settings related to generation of the truth-level event, via Pythia8
        'process' : 'Top_Wqq', # Filepath for a Pythia8 configuration file. Can also correspond to a filename in the `util/pythia/pythia_templates` subdirectory. If no extension is given, `.txt` will be assumed.
        'hadronization' : True, # Pythia8 hadronization flag
        'mpi' : True, # Pythia8 multi-parton interactions flag
        'isr' : True, # Pythia8 initial-state radiation flag
        'fsr' : True, # Pythia8 final-state radiation flag
        'rng' : 1, # Pythia8 RNG seed
        'verbose' : False, # Verbosity flag
        'hepmc_dir': None, # Directory containing the HepMC3 installation. If `None`, will build in a local directory "external/hepmc". Note that we currently use a custom fork of HepMC3, so you probably want to leave this as `None`.
        'hepmc_format': 'root', # Options are 'root' and 'ascii'. The ROOT option provides superior filesize and random-access capability (useful if making samples for pileup overlay), at the cost of being less human-readable.
        'event_filter' : None, # Optional algorithm for filtering out events at generation; for example, require events to have a truth-level jet passing some pT threshold. Events that fail this filter will be replaced, i.e. if you request N events, you will get N events passing this filter (thus one must also be cautious to set it to something reasonable or this may run forever).
    },

    'pileup' : { # settings related to overlay of pileup
        'handler': pu.PileupOverlay( # A class for handling pileup overlay; there are a couple different options. Can set to `None` to skip pileup overlay.
            "/path/to/pileup/hepmc3/*.root", # input files for pileup -- should generally be HepMC3/ROOT format (some classes *require* this format)
            rng_seed=1 # RNG seed for pileup handling
        )
    },

    'simulation' : { # settings related to (fast) detector simulation
        'type' : 'delphes', # what simulation (if any) to use. Currently supported options are `[None, 'delphes']`
        'delphes_card' : "util/delphes/cards/delphes_card_CMS_custom.tcl", # path to the Delphes card to use. If None, will use the ATLAS Delphes card that ships with Delphes.
        'delphes_dir' : None, # Directory containing the Delphes installation. If `None`, will be build in a local directory "external/delphes". Note that our custom CMS card requires a custom fork of Delphes (which will be installed if `None`).
        'delphes_output' : ['EFlowPhoton','EFlowNeutralHadron','EFlowTrack','Electron','Muon','Photon','GenMissingET','MissingET','GenVertex','Vertex'], # Which output objects from Delphes to propagate to the final HDF5 file -- this is also what will be available to the post-processors; other information will be dropped. Note: some details for the vertex-type objects may still need some ironing out.
        'delphes_rng_seed' : 1 # RNG seed for Delphes
    },

    # NOTE: Object names (keys) should not have periods (".") in them. These are used internally to indicate objects' properties ("leaves" in ROOT-speak), and including these in names may break stuff down-the-line (such as in the visualization scripts).
    'reconstruction' : { # settings related to reconstruction (e.g. jet clustering) and n-tupling
        'n_stable' : 200, # max number of stable truth-level particles to save per event (HDF5 doesn't support jagged arrays)
        'n_delphes': [200], # max number of Delphes objects to save per event -- list corresponding to entries in 'delphes_output'. If single value, will be broadcast to the appropriate shape.
        'fastjet_dir' : None, # Directory containing the Fastjet installation. If None, will build in a local directory "external/fastjet". Note that the Fastjet installation you use *must* have the Python bindings set up.
        'n_truth' : 3 + 120, # Maximum number of truth particles to save per event. This should typically be harmonized with the truth-level particle selections set further below. (HDF5 doesn't support jagged arrays)
        'particle_selection' : { # Here, you can specify collections of truth-level particles to save, using various (provided) particle selection algorithms. These typically search for particles matching some PdgID and/or generator status.
            'TruthParticlesTopAndChildren': # this key will be used for naming this collection of particles in the n-tuple
            parsel.MultiSelection( # selection algorithm that consists of a collection of other selection algorithms
                [
                    parsel.FirstSelector(22, 6), # top quark (selects first instance of particle with status code = 22 and PdgID = 6)
                    parsel.FirstSelector(23, 5), # bottom quark
                    parsel.FirstSelector(22,24), # W boson
                    parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(22,24)),n=120) # up to 120 stable daughters of W (stable daughter particles of the first instance of a particle with status code = 22 and pdgID = 24)
                ]
            ),
            'TruthParticlesAntiTopAndChildren': # a second collection of truth particles to save, under this key
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
        'post_processing': [ # What post-processing algorithms to run -- this includes jet clustering! You can queue up multiple separate post-processors. These act on the HDF5 file itself, so you can reference things like truth-particle containers defined above in "particle_selection", or previously-queued post-processors.

            # Cluster large-radius jets, ghost associated to the top and antitop. Enforce some pt and eta cuts.
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08RecoJetsAssociatedTop').PtFilter(25.).EtaFilter(4.).GhostAssociation('TruthParticlesTopAndChildren',    0,mode='filter'), # ghost-associated to top
            jets.JetFinder(['EFlowPhoton','EFlowNeutralHadron','EFlowTrack'],jet_algorithm='anti_kt',radius=0.8,jet_name='AntiKt08RecoJetsAssociatedAntiTop').PtFilter(25.).EtaFilter(4.).GhostAssociation('TruthParticlesAntiTopAndChildren',0,mode='filter'), # ghost-associated to anti-top

            # Cluster small-radius jets, near the ghost-associated jets above -- require them to be within some delta R of the large-radius jets.
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
- `--pythia_config`: A path to the Pythia8 configuration file; this overrides whatever configuration is set in the config file, however this can be a *list* of Pythia8 configuration files in which case `Njobs` jobs will be spun up per Pythia8 configuration. This can be a handy way to simultaneously set up generation of multiple processes. Note however that this may make it difficult to later distinguish between these processes, whereas spinning up separate batches of jobs per process will allow for labeling them via the `signal_flag` setting in the `reconstruction` section of the config file.
- `--shortqueue` : This setting is currently unique to the HTCondor system on the University of Chicago's Analysis Facility, where it allows for accessing a dedicated HTCondor queue for jobs limited to 3 hours. (This code can possibly be extended to similarly access unique features of other HTCondor systems).
- `--n_cpu`: The number of CPU cores to use per job. Note that this code is not (yet) explicitly parallelized, so there is little benefit to requesting multiple cores (aside from that possibly also providing more memory).
- `batch_name`: The name for this batch of condor jobs.
- `mode`: Which "run mode" to use; this is a consequence of HTCondor queues operating in slightly different ways on different clusters. Here are the available options:
  - `1`: The condor jobs will each clone the `HEPData4ML` repository. This requires them to have the ability to clone repositories from GitHub, and any uncommitted local changes will not be present (the clone will use whatever branch of the code you are currently sitting on).
  - `2`: A git clone of the repository will be performed locally, and then this will be shipped to the condor jobs  as a `tar.gz` archive using HTCondor's file transfer protocol.
  - `3`: The code will run "locally", i.e. the condor workers will run the code directly from this repository. This requires the workers to have read access to wherever this repository is located. This mode is useful for example on Brown University's BRUX cluster,  where the HTCondor jobs are not "shipped off" to scratch directories by default.
  - For any other value, the local repository will be packaged into a `tar.gz` archive, and shipped to the condor jobs using HTCondor's file transfer protocol. This is the "typical" run mode.
- `requirements`: A string of requirements for the condor worker nodes.
- `blacklist`: Specific condor worker machines to avoid running on; this can be a handy feature for systems where certain worker nodes don't function well, e.g. some nodes not having consistent access to CVMFS.

**Note** that at the present, the condor workers must have write access to the output directory, as the job explicitly copies output there (instead of using HTCondor's file transfer protocol for output). This is something that may be changed in the future.

### Running with pileup
One of the special features of this package is the ability to include pileup. This is achieved by overlaying a set of "pileup" HepMC3 events with the "main" HepMC3 event before the detector simulation step; it is up to the user to determine how exactly the pileup events are generated. In practice, it is best to generate the pileup events with a dedicated run of this package, invoking `run.py` or `prep_condor.py` with the command line argument
```
--steps "generation"
```
as only the HepMC3 files need to be produced for the pileup. Note that in the config file, the output format should be set to `root` and not `ascii`, as the pileup mixing procedure will rely on random access of the events (which the ASCII format does not support). The ROOT-based format is also more space efficient and generally preferrable (with the caveat being that it is less human-readable than ASCII). Pileup HepMC3 files can of course be supplied externally as well, but those produced with this package have the advantage that they will have metadata attached to them (see further below), which will propagate to the metadata of the dataset they are mixed into; this will greatly help with dataset documentation/reproducibility.

## Parsing and handling output
The final output format of the data generation is an HDF5 file, that contains an n-tuple representing the data. It comprises of a set of n-dimensional `numpy` arrays, which represent properties of the event; the four-momenta of truth-level particles, the four-momenta (and position, charge etc.) of various reconstructed objects from Delphes as well as jet clustering.

HDF5 has the advantage that it is quite widely supported (and relatively popular within the ML community). However unlike ROOT it does not support jagged arrays, which are a relatively natural representation for a lot of HEP data as (for example) the number of particles per event, or number of constituents per jet, is a variable quantity. Thus the n-tuples are constructed using fixed-length arrays, where the array length is (hopefully!) chosen to be sufficiently large as to avoid any edge effects, and zero-padding is employed as needed; each object in the n-tuple also has a column indicating its "multiplicity", i.e. how many of them there are per event.

#### Metadata
In addition to the kinematics and properties of objects in each event, there are also fields corresponding with *metadata*, such as the configuration used to generate the events, with these fields indexing each event with respect to a list of metadata values stored in the HDF5 file's "attributes". When files are concatenated using the Python script `util/tools/concat.py`, merging of the metadata is handled automatically. As it includes the configuration, command line arguments and git hash (among other pieces of information), this metadata in principle allows one to completely reproduce a dataset (aside from the Delphes RNG -- which may be handled in future updates!).

#### Checking dataset contents
One can check the particular contents of a dataset by using the Python script `util/tools/check_file.py`, which will print all the available fields, and their dimensions. This can also be used to print the actual contents of a single event, or check the contents of the metadata, or even print out all the citations (for used software packages and algorithsm) saved in the dataset's metadata in BibTex format:
```
usage: check_file.py [-h] -i INPUTFILE [-ei EVENTINDEX] [-md] [-citations]

options:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
  -ei EVENTINDEX, --eventIndex EVENTINDEX
  -md, --metaData       Check the metadata values.
  -citations, --citations
                        Print list of all citations for algorithms used.
```

