# This sets up the environment with a Python3 installation with most dependencies (numpy, h5py, pandas, uproot, PyROOT). Importantly, it also has the Pythia8 Python bindings.
source /cvmfs/sft.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc12-opt/setup.sh # see https://sft.its.cern.ch/jira/projects/GENSER/issues/GENSER-440

# We will still need to use pip to install PyHepMC locally, this is not included in the above.
pip install --user pyhepmc
