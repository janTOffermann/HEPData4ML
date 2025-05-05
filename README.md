# HEPData4ML

This package is meant to provide a relatively easy way to generate Monte Carlo datasets for high-energy physics (HEP) studies, such as (but not limited to) development of machine learning methods.

## Overview

Specifically, this software package is designed for generating datasets for studying jets, although it may be extended to other uses. It uses the [Pythia8](https://pythia.org) Monte Carlo generator to simulate proton-proton collisions (or *events*), and applies a set of user-defined filters to determine what information to save from these events. It then (optionally) pipes this through fast detector simulation via [DELPHES](https://cp3.irmp.ucl.ac.be/projects/delphes), and performs jet clustering on the output using [fastjet](http://fastjet.fr). It then selects one jet from each event, and saves this jet along with its constituents, as well as any user-selected "truth-level" particles from the event.

In effect, this package is a Pythia (+Delphes) wrapper that produces n-tuples in HDF5 format. Its design makes it most useful for jet physics -- as there is a jet clustering step in the workflow described above, and each output entry corresponds not to a full event but a single jet -- but with some advanced usage the jet clustering and selection steps can be omitted, allowing for general production of Pythia8 events (possibly with detector simulation).

## Requirements/setup

The code is primarily written in Python3, along with some C++/ROOT code for fast calculations involving four-vectors. It is compatible with both Linux and macOS operating systems -- for the latter, only Intel-based Macs are currently supported (since fastjet does not currently support compilation on ARM-based Macs).

HEPData4ML makes use of a number of different software packages, including:

- [ROOT](https://root.cern.ch)
- [numpy](https://numpy.org)
- [h5py](https://www.h5py.org)
- [pyhepmc](https://github.com/scikit-hep/pyhepmc)
- [Pythia8](https://pythia.org)
- [DELPHES](https://cp3.irmp.ucl.ac.be/projects/delphes)
- [fastjet](http://fastjet.fr)

The easiest way to prepare the software environment is via the [conda package manager](https://docs.conda.io/en/latest/). A YAML file detailing the necessary conda environment is provided in `setup/conda/`. This includes a binary distribution of ROOT as provided by conda-forge. The DELPHES and fastjet dependencies will automatically be installed locally (i.e. *within the repo*, to make removal easy) by the package during runtime, if not found.

## Configuring the data generation
The majority of the settings for event generation are set within `util/config.py`, within a dictionary called `config`. Here is an overview of its entries:

- `proc`: A physics process to run. The string must correspond with the name of a Pythia8 configuration file in `util/pythia_templates/`, wherein the Pythia8 processes can be specified. *Note that a few Pythia8 configurations are separately handled by other entries in* `config`.
- `hadronization`: Whether or not to enable hadronization in Pythia8. (the `HadronLevel:all` Pythia8 option)
- `mpi`: Whether or not to enable multi-parton interactions in Pythia8. (the `PartonLevel:MPI` Pythia8 option)
- `isr`: Whether or not to enable initial state radiation in Pythia8. (the `PartonLevel:ISR` Pythia8 option)
- `fsr`: Whether or not to enable final state radiation in Pythia8. (the `PartonLevel:FSR` Pythia8 option)
- `delphes`: Whether or not to enable DELPHES (fast detector simulation). If enabled, DELPHES fast detector simulation (using the ATLAS detector card) is performed, and its output `Towers` are passed to FastJet for clustering. If disabled, we simply pass the final-state (stable) particles from Pythia8 output to FastJet.
- `rng`: The random number generator seed used by Pythia8.
- `jet_radius`: The jet radius to be used by FastJet (using the anti-kt algorithm).
- `jet_pt_min`: A minimum pT cut selection to apply to jets.
- `jet_max_eta`: A maximum abs(eta) cut selection to apply to jets.
- `jet_n_par`: The maximum number of jet constituents to save for each jet (the list of constituents will be zero-padded).
- `n_truth`: The number of truth-level particles to save for each jet. *This will be further explained below.*
- `truth_selection`: An algorithm for selecting which truth-level particles to save for each jet.
- `final_state_selection`: An algorithm for selecting which final-state (stable) particles to pass on to DELPHES or FastJet.
- `event_selection`: An (optional) algorithm for further filtering the final-state particles above.
- `jet_selection`: An algorithm for selecting which jets to choose from each Pythia8 event.
- `signal_flag`: An integer to mark these events -- this is useful if producing multiple datasets that will be combined, so that one can label different datasets as signal or background.
- `split_seed`: The RNG seed used when splitting the dataset into training, validation and testing samples.

The selectors `config` (`truth_selection`, `final_state_selection`, `event_selection`, `jet_selection`) correspond with how final-state particles, truth-level particles and jets are chosen from each Pythia8 event:

### Truth-level selection

The `truth_selection` algorithm will determine what particles we save as our "truth-level" information. This will not be passed to jet clustering, but will be present in the final dataset output. These can be any particles from the Pythia8 event listing, and is typically useful for storing classification/regression targets. For example, if producing top quark jets, one may save the truth-level top quark, as well as the bottom quark and W boson it decays to.

### Final-state selection

The `final_state_selection` algorithm will determine what particles we save as our final state. If using DELPHES, this is what will be passed to DELPHES as input[^1]. If not, this is what will be directly passed to jet clustering. For example, we may just want to save all stable particles output by Pythia.

[^1]: Of course, DELPHES will by default only actually interact with certain particles by default (e.g. the detector does not interact with neutrinos). Nonetheless, this may prove a useful handle for deciding what particles propagate forward in the simulation process.

### Event selection

The `event_selection` algorithm is perhaps better described as a filter -- this lets us further filter the output of `final_state_selection`, which may be handy if we know there are certain particles there that we don't actually care about. For example, if producing top quark jets via ttbar processes, we will have two top quark jets per event (one from the top quark, one from the top anti-quark). If we're selecting the truth-level top we are probably not interested in the jet produced by the anti-top, so if `final_state_selection` has picked up all stable particles we may want to throw away those that are very far from our truth-level top (which are possibly just produced by the anti-top, and which won't go into our top quark jet anyway if they are sufficiently far away). The `event_selection` can be set to `None` to avoid this additional filtering.

### Jet selection

The `jet_selection` algorithm is used to determine which jet to select for each event -- as a current design limitation, we only select one jet per event (and discard the others). *In principle, one can set this to `None` in order to skip jet clustering, in which case the dataset will just contain all the final-state particles instead of jet constituents*.

## Running the code

With the configuration set, the code can then be run by invoking `run.py` as follows:

```
usage: run.py [-h] -n NEVENTS -p PTBINS [PTBINS ...] [-o OUTFILE] [-O OUTDIR]
              [-g GENERATION] [-s SEP_TRUTH] [-ns N_SEP_TRUTH]
              [-d DIAGNOSTIC_PLOTS] [-v VERBOSE] [-h5 HDF5] [-f FORCE]
              [-c COMPRESS] [-cd CLEAN_DELPHES]

optional arguments:
  -h, --help            show this help message and exit
  -n NEVENTS, --nevents NEVENTS
                        Number of events per pt bin.
  -p PTBINS [PTBINS ...], --ptbins PTBINS [PTBINS ...]
                        Transverse momentum bin edges.
  -o OUTFILE, --outfile OUTFILE
                        Output HDF5 file name.
  -O OUTDIR, --outdir OUTDIR
                        Output directory.
  -g GENERATION, --generation GENERATION
                        Whether or not to do event generation.
  -s SEP_TRUTH, --sep_truth SEP_TRUTH
                        Whether or not to store truth-level particles in
                        separate arrays.
  -ns N_SEP_TRUTH, --n_sep_truth N_SEP_TRUTH
                        How many truth particles to save in separate arrays --
                        will save the first n as given by the truth selection.
  -d DIAGNOSTIC_PLOTS, --diagnostic_plots DIAGNOSTIC_PLOTS
                        Whether or not to make diagnostic plots.
  -v VERBOSE, --verbose VERBOSE
                        Verbosity.
  -h5 HDF5, --hdf5 HDF5
                        Whether or not to produce final HDF5 files. If false,
                        stops after HepMC or Delphes/ROOT file production.
  -f FORCE, --force FORCE
                        Whether or not to force generation -- if true, will
                        possibly overwrite existing HepMC files in output
                        directory.
  -c COMPRESS, --compress COMPRESS
                        Whether or not to compress HepMC files.
  -cd CLEAN_DELPHES, --clean-delphes CLEAN_DELPHES
                        Whether or not to clean up DELPHES/ROOT files.
```
Here are a few notes on the arguments passed to `run.py`:

- As noted above, `NEVENTS` indicates the number of events to generate *per pT bin*. Note that generating an event does not necessarily mean that it will have any jets that pass the pT and eta cuts, or the `jet_selection` algorithm, so the final number of jets per bin may be smaller than `NEVENTS`.
- The values of `PTBINS` correspond with the "pT hat" parameter used by Pythia8 -- this gives the pT of the leading "parton final-state" particle (i.e. an outgoing particle from the selected tree-level process, *before* any decay or hadronization). These values are given in GeV.
- The `GENERATION` flag can be used to skip event generation, which can be done *if* the HEPMC files (which are generated by this routine) are present. This can be used to quickly remake the final HDF5 files, as the event generation (which produces the HEPMC files) is typically the most time-consuming step.