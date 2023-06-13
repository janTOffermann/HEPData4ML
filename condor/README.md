##Condor usage

This directory contains resources for running our generation using the [HTCondor](https://htcondor.org) batch system.

### prep_condor.py
This script will prepare a bunch of jobs for condor. Specifically it will package up the generation software, and prepare a condor submission file together with the script that the condor job will run, as well as a plaintext file that provides the arguments for each job. Since it's difficult to cover all possible use cases of the parallelization without making this really convoluted, this design allows the user to then modify that plaintext file however they would like, e.g. adjust the pT bins of the various jobs, or make them run with different configurations. Each set of jobs (a condor "cluster" of jobs) will use a single configuration file (config.py) that provides much of the configuration -- everything that is not explicitly set by command line arguments for the `run.py` script.