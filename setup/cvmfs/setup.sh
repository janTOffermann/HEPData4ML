Help()
{
  # Display Help
  echo "Add description of the script functions here."
  echo
  echo "Syntax: scriptTemplate [-h|l]"
  echo "options:"
  echo "h     Print this Help."
  echo "l     Force CentOs7 instead of EL9 setup."
  echo
}

# ============================
#        Main program
# ============================

lcg="LCG_105"
setupOption="EL9"
build="x86_64-el9-gcc13-opt"

while getopts ":hl" option; do # getopts is kind of terrible, but this will do for now...
  case $option in
    h) # display Help
      Help
      return;;
    l)
      setupOption="CentOS7"
      ;;
  esac
done

# # The following couple of lines are the equivalent of the "setupATLAS" command that is commonly used in ATLAS for environment setup.
# if [ -z $ATLAS_LOCAL_ROOT_BASE ]; then
#   export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
# fi
# source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh -q

echo "Checking OS version. (see setup.sh script)"
os_version=$(uname -r)
echo "OS version: ${os_version}"
# if [[ $os_version == *"el9"* ]]; then
#   setupOption="EL9"
# else
#   setupOption="Centos7"
# fi

echo "Setup for ${setupOption}"

# # Try to automatically check OS version, and adjust setupOption. (could ultimately remove the user setting.)
# TODO: This breaks condor jobs when EL9 is detected:
#
# /cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc13-opt/lib64/R/bin/exec/R: error while loading shared libraries: libreadline.so.8: cannot open shared object file: No such file or directory
# /cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc13-opt/lib64/R/bin/exec/R: error while loading shared libraries: libreadline.so.8: cannot open shared object file: No such file or directory
# /cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc13-opt/lib64/R/bin/exec/R: error while loading shared libraries: libreadline.so.8: cannot open shared object file: No such file or directory
# /cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc13-opt/lib64/R/bin/exec/R: error while loading shared libraries: libreadline.so.8: cannot open shared object file: No such file or directory
# /cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc13-opt/lib64/R/bin/exec/R: error while loading shared libraries: libreadline.so.8: cannot open shared object file: No such file or directory
# gnuplot: error while loading shared libraries: libreadline.so.8: cannot open shared object file: No such file or directory
# python: error while loading shared libraries: libcrypt.so.2: cannot open shared object file: No such file or directory
#

# os_version=$(uname -r)

# if [[ $os_version =~ "el7" ]]; then
#   setupOption="Centos7"
#   echo "Detected EL7, adjusting setup accordingly."
# fi

# if [[ $os_version =~ "el9" ]]; then
#   setupOption="EL9"
#   echo "Detected EL9, adjusting setup accordingly."
# fi

# --------------------------------------------------

if [ "$setupOption" == "CentOS7" ]; then
  build="x86_64-centos7-gcc11-opt"
fi
echo "Using ${lcg} with build ${build}."

# For ATLAS, was using "lsetup" but this should be more portable.
# (It also is consistent with software setup in CMS Data Analysis School tutorials)
source /cvmfs/sft.cern.ch/lcg/views/${lcg}/${build}/setup.sh

#echo "-----------------------"
#echo "Here is a potentially useful option for the configuration file in config/config.py:"
#echo "'delphes_dir' : '/cvmfs/sft.cern.ch/lcg/releases/delphes/3.5.1pre05-775ca/x86_64-centos7-gcc11-opt'"
#echo "For fastjet, you will need Python bindings and thus a local build -- preferably done interactively, so that you can later point condor jobs to that and they don't all have to build fastjet again."
