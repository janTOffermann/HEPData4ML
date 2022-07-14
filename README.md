# HEPData4ML

This package is meant to provide a relatively easy way to generate Monte Carlo datasets for high-energy physics (HEP) studies, such as (but not limited to) development of machine learning methods. Specifically, this software package is designed for generating datasets for studying jets, although it may be extended to other uses.

## Requirements/setup

The code is primarily written in Python3, along with some C++/ROOT code for fast calculations involving four-vectors. It is compatible with both Linux and macOS operating systems.

HEPData4ML makes use of a number of different software packages, including:

- [ROOT](https://root.cern.ch)
- [numpy](https://numpy.org)
- [h5py](https://www.h5py.org)
- [pyhepmc-ng](https://github.com/scikit-hep/pyhepmc)
- [numpythia](https://github.com/scikit-hep/numpythia)
    - [Pythia8](https://pythia.org) (used implicitly via numpythia)
- [DELPHES](https://cp3.irmp.ucl.ac.be/projects/delphes)
- [fastjet](http://fastjet.fr)

The easiest way to prepare the software environment is via the [conda package manager](https://docs.conda.io/en/latest/). A YAML file detailing the necessary conda environment is provided in `setup/conda/`. This includes a binary distribution of ROOT as provided by conda-forge. The DELPHES and fastjet dependencies will automatically be installed locally (i.e. *within the repo*, to make removal easy) by the package during runtime, if not found.

## Configuring the data generation
The majority of the settings for event generation are set within `util/config.py`, within a dictionary called `config`. Here is an overview of its entries:

- `proc`: A physics process to run. The string must correspond with the name of a Pythia8 configuration file in `util/pythia_templates/`, wherein the Pythia8 processes can be specified. *Note that a few Pythia8 configurations are separately handled by other entries in* `config`.
- `hadronization`: Whether or not to enable hadronization in Pythia8.
- `mpi`: Whether or not to enable multi-parton interactions in Pythia8.
- `isr`: Whether or not to enable initial state radiation in Pythia8.
- `fsr`: Whether or not to enable final state radiation in Pythia8.
- `delphes`: Whether or not to enable DELPHES (fast detector simulation). If off, we don't perform any detector simulation (we just use Pythia8 output).
- `rng`: The random number generator seed used by Pythia8.
- `jet_radius`: The jet radius to be used by FastJet (using the anti-kt algorithm).
- `jet_pt_min`: A minimum pT cut selection to apply to jets.
- `jet_max_eta`: A maximum abs(eta) cut selection to apply to jets.
- `jet_n_par`: The maximum number of jet constituents to save for each jet (the list of constituents will be zero-padded).
- `n_truth`: The number of truth-level particles to save for each jet. *This will be further explained below.*
- `truth_selection`: An algorithm for selecting which truth-level particles to save for each jet.
- `jet_selection`: An algorithm for selecting which jets to choose from each Pythia8 event.

The last three entries in `config` (`n_truth`, `truth_selection`, `jet_selection`) correspond with how jets and truth-level particles are chosen from each Pythia8 event: An event may have multiple jets, but we (currently) only select one jet from each event, with the selection method determined by whatever function is assigned to `jet_selection`. Similarly, we may choose to save some truth-level particles from the Pythia8 event listing, with the number of particles given by `n_truth` and the way they are selected given by `truth_selection`. (With the current design, we expect to save the same number of truth-level particles for each event.)

The `truth_selection` and `jet_selection` algorithms are defined in `util/truth_selection.py` and `util/jet_selection.py` accordingly.

## Running the code

With the configuration set, the code can then be run by invoking `run.py` as follows:

```
usage: run.py [-h] -n NEVENTS -p PTBINS [PTBINS ...] [-o OUTFILE] [-O OUTDIR]
              [-g GENERATION] [-s SEP_TRUTH]

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
```
Here are a few notes on the arguments passed to `run.py`:

- As noted above, `NEVENTS` indicates the number of events to generate *per pT bin*. Note that generating an event does not necessarily mean that it will have any jets that pass the pT and eta cuts, or the `jet_selection` algorithm, so the final number of jets per bin may be smaller than `NEVENTS`.
- The values of `PTBINS` correspond with the "pT hat" parameter used by Pythia8 -- this gives the pT of the leading "parton final-state" particle (i.e. an outgoing particle from the selected tree-level process, *before* any decay or hadronization). These values are given in GeV.
- The `GENERATION` flag can be used to skip event generation, which can be done *if* the HEPMC files (which are generated by this routine) are present. This can be used to quickly remake the final HDF5 files, as the event generation (which produces the HEPMC files) is typically the most time-consuming step.