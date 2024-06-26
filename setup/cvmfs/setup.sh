# This sets up the environment with a Python3 installation with most dependencies (numpy, h5py, pandas, uproot, PyROOT). Importantly, it also has the Pythia8 Python bindings.
# This is a "typical" ATLAS software setup. Commented out further below, there is a possibly more fragile (but not ATLAS-specific) way to get most of the software,
# but it does not include fastjet (which this package will install locally if needed).

# The following couple of lines are the equivalent of the "setupATLAS" command that is commonly used in ATLAS for environment setup.
if [ -z $ATLAS_LOCAL_ROOT_BASE ]; then
        export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
fi
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh -q

# Now we set up a "view" of an LCG release. The LCG_103 contains Pythia8 bindings, earlier ones may not.
# See: https://sft.its.cern.ch/jira/projects/GENSER/issues/GENSER-440
lsetup "views LCG_105 x86_64-el9-gcc13-opt" # may technically need to adjust the last part based on the arch of where you run this

#echo "-----------------------"
#echo "Here is a potentially useful option for the configuration file in config/config.py:"
#echo "'delphes_dir' : '/cvmfs/sft.cern.ch/lcg/releases/delphes/3.5.1pre05-775ca/x86_64-centos7-gcc11-opt'"
#echo "For fastjet, you will need Python bindings and thus a local build -- preferably done interactively, so that you can later point condor jobs to that and they don't all have to build fastjet again."
