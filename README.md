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

The code is primarily written in Python3, along with some C++/ROOT code (currently only leveraged for the Johns Hopkins top tagger implementation from fastjet).

HEPData4ML makes use of a number of different software packages, including:

- [ROOT](https://root.cern.ch)
- [uproot v4](https://uproot.readthedocs.io/en/latest/index.html)
  - [awkward v1](https://awkward-array.org/doc/main/)
- [numpy](https://numpy.org)
- [h5py](https://www.h5py.org)
- [pyhepmc](https://github.com/scikit-hep/pyhepmc)
- [Pythia8](https://pythia.org)
- [DELPHES](https://cp3.irmp.ucl.ac.be/projects/delphes)
- [fastjet](http://fastjet.fr)

A YAML file detailing a `conda` environment with most of these dependencies is available in `setup/conda/`, including a binary distribution of ROOT.
For usage on CVMFS systems, the setup script in `setup/cvmfs/` is recommended instead.
The DELPHES and fastjet dependencies will automatically be installed locally (i.e. within this repository) by the package during runtime, if they are not found; their locations should be provided in the aforementioned Python configuration file.

Note that the Pythia8 and fastjet installations that are used must have the Python interface built.

## Configuring the data generation

**TODO** The rest of this readme needs to be rewritten; the code on the `devel_full_event` branch is cleaner than on `main` but differs appreciably in how it is run (there are fewer command-line arguments, and the configuration file is a little more organized).