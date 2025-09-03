#!/bin/bash
# This is the shell script that will be run as the condor job.

input_file=$1
proc_number=$2
git_option=$3
git_branch=$4

local_mode=0

# Set up the code. This may involve shipping in a payload, or running `git clone` here.
gitdir=HEPData4ML # TODO: payload curently set to use this name, is this OK or too much hardcoding?
if [[ "${git_option}" == "1" ]]; then
  # run git clone here
  echo "Cloning code from GitHub."
  git clone -b ${git_branch} git@github.com:janTOffermann/HEPData4ML.git ${gitdir}
elif [[ -d "${git_option}" ]]; then
  gitdir=$git_option # the $git_option variable is actually being used to give a path to existing HepData4ML installation
  echo "Running from ${gitdir} ."
  local_mode=1 # need to be a bit careful that we don't delete useful files!
else
  # assume the payload has been shipped in, as payload.tar.gz
  payload=payload.tar.gz
  tar -xzf $payload
  rm $payload
  echo "Unpacked code from payload ${payload} ."
fi

# Run the setup script.
source ${gitdir}/setup/cvmfs/setup.sh

# # Set the number of threads (for OpenBLAS). Might be necessary in order to deal with memory limits.
# export OPENBLAS_NUM_THREADS=${openblas_max_thread}
# export GOTO_NUM_THREADS=${openblas_max_thread}
# export OMP_NUM_THREADS=${openblas_max_thread}

# ===========================================================
echo "Invoking util/tools/compute_pileup_pt.py..."
python ${gitdir}/util/tools/compute_pileup_pt.py \
  --inputFiles $input_file

if [[ "${local_mode}" == "0" ]]; then
  rm -rf $gitdir
  if [ -d external/fastjet ]; then rm -r external/fastjet; fi
  if [ -d external/delphes ]; then rm -r external/delphes; fi
  rm *.py
fi
