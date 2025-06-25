#!/bin/bash
# This is the shell script that will be run as the condor job.
# Outline of the arguments:
# $1 Number of events per pT bin.
# $2 pT bins (list of bin edges)
# $3 which steps to run
# $3 RNG seed for generation. (can be used to overwrite the builtin config file)
# $4 whether or not to split final HDF5 file into train/validation/test files. Only relevant if making the HDF5 file.
# $5 Pythia config (can be used to overwrite the builtin config file)
# $6 Event index offset.
# $7 Configuration file (Python).
# $8 Output directory (for the condor job).
# $9 Process number (for naming the output).
# $10 OpenBLAS max thread count (for multithreading).
# $11 Git option. Determines if we do a git clone here, or if the code has been shipped in as a tarball.
# $12 Git branch.
#

nevents_per_bin=$1
pt_bins=$2
steps=$3
rng_seed=$4
do_split=$5
pythia_config=$6
event_idx_offset=$7
config_file=$8
outdir=$9
proc_number=${10}
openblas_max_thread=${11}
git_option=${12}
git_branch=${13}

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

delete_delphes=0
training_faction=0.6
validation_fraction=0.2

# ===========================================================
echo "Invoking run.py..."
python ${gitdir}/run.py \
  -n ${nevents_per_bin} \
  -p ${pt_bins} \
  -steps ${steps} \
  -O ${outdir_local} \
  -o ${output_filename} \
  -rng ${rng_seed} \
  -pb 1 \
  --delete_stats 0 \
  --split 0 \
  -del_delphes $delete_delphes \
  -pc ${pythia_config} \
  -df 1 \
  --index_offset ${event_idx_offset} \
  --config ${config_file}

copy_script=${gitdir}/util/condor/copy_output.py

if [ "${do_split}" == "1" ]; then
  echo "Splitting output."

  python ${gitdir}/util/tools/split.py \
    -i $output_file \
    -o ${outdir_local} \
    -f1 $training_faction \
    -f2 $validation_fraction \
    -s 1 \
    -c 9
  rm $output_file
  python $copy_script -i ${outdir_local}/train.h5 -e "h5" -o ${outdir} -n ${proc_number}
  python $copy_script -i ${outdir_local}/test.h5  -e "h5" -o ${outdir} -n ${proc_number}
  python $copy_script -i ${outdir_local}/valid.h5 -e "h5" -o ${outdir} -n ${proc_number}
  rm ${outdir_local}/train.h5 ${outdir_local}/test.h5 ${outdir_local}/valid.h5

else
  python $copy_script -i $output_file -e "h5" -o ${outdir} -n ${proc_number}
  rm $output_file
fi

# Compress the output and extract it.
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
