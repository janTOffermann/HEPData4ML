#!/bin/bash
# This is the shell script that will be run as the condor job.
# Outline of the arguments:
# $1 Number of events per pT bin.
# $2 pT bins (list of bin edges)
# $3 which steps to run
# $4 RNG seed for generation. (can be used to overwrite the builtin config file)
# $5 Pythia config (can be used to overwrite the builtin config file)
# $6 Event index offset.
# $7 Job number (TODO: redundant with Process number).
# $8 Total number of jobs.

# $9 Configuration file (Python).
# $10 Output directory (for the condor job).
# $11 Process number (for naming the output).
# $12 OpenBLAS max thread count (for multithreading).
# $13 Git option. Determines if we do a git clone here, or if the code has been shipped in as a tarball.
# $14 Git branch.
#

nevents_per_bin=$1
pt_bins=$2
steps=$3
rng_seed=$4
pythia_config=$5
event_idx_offset=$6
job_number=$7
njobs_total=$8
config_file=$9
outdir=${10}
proc_number=${11}
openblas_max_thread=${12}
git_option=${13}
git_branch=${14}

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

# Set the number of threads (for OpenBLAS). Might be necessary in order to deal with memory limits.
export OPENBLAS_NUM_THREADS=${openblas_max_thread}
export GOTO_NUM_THREADS=${openblas_max_thread}
export OMP_NUM_THREADS=${openblas_max_thread}

# Move the config.py file into the config directory. It has been shipped as an input file separate of the payload.
# TODO: Could be an issue for local_mode=1 -- why wasn't I previously passing $config_file as an arg below? Maybe will rediscover some old bug.
# mv $config_file ${gitdir}/config/config.py

outdir_local="output_${proc_number}"
output_filename=events.h5
output_file="${outdir_local}/${output_filename}"

# ===========================================================
echo "Invoking run.py..."
python ${gitdir}/run.py \
  -n ${nevents_per_bin} \
  --ptbins="${pt_bins}" \
  -steps ${steps} \
  -O ${outdir_local} \
  -o ${output_filename} \
  -rng ${rng_seed} \
  -pb 1 \
  -pc ${pythia_config} \
  --index_offset ${event_idx_offset} \
  --config ${config_file} \
  --condor \
  --condor_job_number ${job_number} \
  --n_condor_jobs ${njobs_total}

copy_script=${gitdir}/util/condor/copy_output.py

# NOTE: The output_file doesn't necessarily exist; if we to the generation step only
#       we won't have it. (We will have some other output, that's in the outdir_local
#       and that we'll capture further below.

if test -f "${output_file}"; then
  python $copy_script -i $output_file -e "h5" -o ${outdir} -n ${proc_number}
  rm $output_file
fi

# Compress the full output and extract it.
outname="output.tar.gz"
tar -czf ${outname} ${outdir_local}

# Ship the output tarball.
python $copy_script -i ${outname} -e "tar.gz" -o ${outdir} -n ${proc_number}

# Cleanup. Not strictly necessary.
rm -r ${outdir_local}
rm ${outname}

if [[ "${local_mode}" == "0" ]]; then
  rm -rf $gitdir
  if [ -d external/fastjet ]; then rm -r external/fastjet; fi
  if [ -d external/delphes ]; then rm -r external/delphes; fi
  rm *.py
fi
