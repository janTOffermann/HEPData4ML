#!/bin/bash
# This is the shell script that will be run as the condor job.
# Outline of the arguments:
# $1 Number of events per pT bin.
# $2 pT bins (list of bin edges)
# $3 whether or not to save separate truth particle keys (0 or 1)
# $4 number of separate truth particles to save (i.e. save the first n from the list)
# $5 whether or not to do jet clustering and make HDF5 files (HepMC -> HDF5). If not, stops after making the HepMC files.
# $6 RNG seed for generation. (can be used to overwrite the builtin config file)
# $7 whether or not to split final HDF5 file into train/validation/test files. Only relevant if making the HDF5 file.
# $8 Pythia config (can be used to overwrite the builtin config file)
# $9 Event index offset.
# $10 Configuration file (Python).
# $11 Output directory (for the condor job).
# $12 Process number (for naming the output).
# $13 Whether or not to run Delphes. If yes (>0), will run both truth-level and Delphes, and match them so that they have the same events.
# $14 OpenBLAS max thread count (for multithreading).
# $15 Git option. Determines if we do a git clone here, or if the code has been shipped in as a tarball.
# $16 Git branch.

nevents_per_bin=$1
pt_bins=$2
separate_truth_flag=$3
separate_truth_number=$4
do_h5=$5
rng_seed=$6
do_split=$7
pythia_config=$8
event_idx_offset=$9
config_file=${10}
outdir=${11}
proc_number=${12}
do_delphes=${13}
openblas_max_thread=${14}
git_option=${15}
git_branch=${16}

# Set up the code. This may involve shipping in a payload, or running `git clone` here.
gitdir=HEPData4ML # TODO: payload curently set to use this name, is this OK or too much hardcoding?
if [ "${git_option}" == "0" ]; then
  # run git clone here
  git clone -b ${git_branch} git@github.com:janTOffermann/HEPData4ML.git ${gitdir}
else
  # assume the payload has been shipped in, as payload.tar.gz
  payload=payload.tar.gz
  tar -xzf ${payload}
  rm ${payload}
fi

# Run the setup script.
source ${gitdir}/setup/cvmfs/setup.sh

# Set the number of threads (for OpenBLAS). Might be necessary in order to deal with memory limits.
export OPENBLAS_NUM_THREADS=${openblas_max_thread}
export GOTO_NUM_THREADS=${openblas_max_thread}
export OMP_NUM_THREADS=${openblas_max_thread}

# Move the config.py file into the config directory. It has been shipped as an input file separate of the payload.
mv $config_file ${gitdir}/config/config.py

outdir_local="output_${proc_number}"
truth_file_short=events.h5
delphes_file_short=events_delphes.h5
truth_file="${outdir_local}/${truth_file_short}"
delphes_file="${outdir_local}/${delphes_file_short}"

# ===========================================================
# If Delphes was requested, we run the "standard" workflow:
#  1) Run with Delphes.
#  2) Run without Delphes (i.e. truth-level).
#  3) Filter so Delphes & truth have only the same events.
#     - this is because any jet-level selections placed
#       on Delphes or truth-level jets may result in
#       different events getting discarded from each run.
# ===========================================================
if [ "$do_delphes" -gt "0" ]; then
  # Run the generation. First with Delphes, then again without.
  # (for the second run, use the HepMC files that have already been generated)
  echo "Generating events, applying DELPHES fast detector sim, clustering jets."
  python ${gitdir}/run.py \
    -n ${nevents_per_bin} \
    -p_s ${pt_bins} \
    -O ${outdir_local} \
    -o ${delphes_file_short} \
    -s ${separate_truth_flag} \
    -ns ${separate_truth_number} \
    -h5 ${do_h5} \
    -rng ${rng_seed} \
    -pb 1 \
    --delete_stats 0 \
    --split 0 \
    -del_delphes 1 \
    -pc ${pythia_config} \
    -df 1 \
    --index_offset ${event_idx_offset} \
    --delphes 1
    # --config ${config_file} \

  # Run a 2nd time, now without Delphes output. Produces truth-level "events.h5"
  echo "Now clustering truth-level jets from previously-generated events."
  python ${gitdir}/run.py \
    -n ${nevents_per_bin} \
    -p_s ${pt_bins} \
    -O ${outdir_local} \
    -o ${truth_file_short} \
    -s ${separate_truth_flag} \
    -ns ${separate_truth_number} \
    -h5 ${do_h5} \
    -rng ${rng_seed} \
    -pb 0 \
    --delete_stats 0 \
    --split 0 \
    -pc ${pythia_config} \
    -df 1 \
    --index_offset ${event_idx_offset} \
    --delphes 0
    # --config ${config_file} \

  # Slim down to only events shared between the truth-level and
  # Delphes files (events that have passed jet cuts for both the truth-level
  # and Delphes jets).
  # TODO: Make this optional?
  truth_file_matched="${outdir_local}/events_matched.h5"
  delphes_file_matched="${outdir_local}/events_delphes_matched.h5"
  python ${gitdir}/util/tools/match.py \
    -i1 $truth_file \
    -i2 $delphes_file \
    -o1 $truth_file_matched \
    -o2 $delphes_file_matched \
    -k event_idx

  mv $truth_file_matched $truth_file
  mv $delphes_file_matched $delphes_file

  # Optionally split files into train,test and validation sets.
  # TODO: Split is currently hard-coded, make this configurable?
  if [ "${do_split}" == "1" ]; then
    echo "Splitting files."
    python ${gitdir}/util/tools/split.py \
      -i $delphes_file \
      -o ${outdir_local} \
      -f1 0.6 \
      -f2 0.2 \
      -s 1 \
      -c 9

    rm $delphes_file
    # Prepend "delphes" to help with possible sorting of outputs (delphes vs. truth),
    # useful when there are many files.
    mv ${outdir_local}/train.h5 ${outdir_local}/delphes_train.h5
    mv ${outdir_local}/test.h5 ${outdir_local}/delphes_test.h5
    mv ${outdir_local}/valid.h5 ${outdir_local}/delphes_valid.h5

    python ${gitdir}/util/tools/split.py \
      -i $truth_file \
      -o ${outdir_local} \
      -f1 0.6 \
      -f2 0.2 \
      -s 1 \
      -c 9

    rm $truth_file
  fi

  # Ship the output HDF5 file(s).
  # We will also ship the other files that are produced, but will bundle them up.
  # (we keep these separate since in practice they are the most useful to directly access)
  if [ "${do_split}" == "1" ]; then
    python copy_output.py -i ${outdir_local}/train.h5 -e "h5" -o ${outdir} -n ${proc_number}
    python copy_output.py -i ${outdir_local}/test.h5  -e "h5" -o ${outdir} -n ${proc_number}
    python copy_output.py -i ${outdir_local}/valid.h5 -e "h5" -o ${outdir} -n ${proc_number}
    rm ${outdir_local}/train.h5 ${outdir_local}/test.h5 ${outdir_local}/valid.h5

    python copy_output.py -i ${outdir_local}/delphes_train.h5 -e "h5" -o ${outdir} -n ${proc_number}
    python copy_output.py -i ${outdir_local}/delphes_test.h5  -e "h5" -o ${outdir} -n ${proc_number}
    python copy_output.py -i ${outdir_local}/delphes_valid.h5 -e "h5" -o ${outdir} -n ${proc_number}
    rm ${outdir_local}/delphes_train.h5 ${outdir_local}/delphes_test.h5 ${outdir_local}/delphes_valid.h5

  else
    python copy_output.py -i $truth_file -e "h5" -o ${outdir} -n ${proc_number}
    python copy_output.py -i $delphes_file -e "h5" -o ${outdir} -n ${proc_number}

    rm $truth_file
    rm $delphes_file
  fi
# ===========================================================
# If Delphes NOT requested, just run the truth-level stuff.
# This is simpler, because there's no need to do any kind
# of event matching between two runs.
# ===========================================================
else
  echo "Generating events (truth-level only, no DELPHES)."
  python ${gitdir}/run.py \
    -n ${nevents_per_bin} \
    -p_s ${pt_bins} \
    -O ${outdir_local} \
    -o ${truth_file_short} \
    -s ${separate_truth_flag} \
    -ns ${separate_truth_number} \
    -h5 ${do_h5} \
    -rng ${rng_seed} \
    -pb 1 \
    --delete_stats 0 \
    --split 0 \
    -del_delphes 1 \
    -pc ${pythia_config} \
    -df 1 \
    --index_offset ${event_idx_offset} \
    --delphes 0
    # --config ${config_file} \

  if [ "${do_split}" == "1" ]; then
    echo "Splitting output."

    python ${gitdir}/util/tools/split.py \
      -i $truth_file \
      -o ${outdir_local} \
      -f1 0.6 \
      -f2 0.2 \
      -s 1 \
      -c 9
    rm $truth_file
    python copy_output.py -i ${outdir_local}/train.h5 -e "h5" -o ${outdir} -n ${proc_number}
    python copy_output.py -i ${outdir_local}/test.h5  -e "h5" -o ${outdir} -n ${proc_number}
    python copy_output.py -i ${outdir_local}/valid.h5 -e "h5" -o ${outdir} -n ${proc_number}
    rm ${outdir_local}/train.h5 ${outdir_local}/test.h5 ${outdir_local}/valid.h5

  else
    python copy_output.py -i $truth_file -e "h5" -o ${outdir} -n ${proc_number}
    rm $truth_file
  fi
fi

# Compress the output and extract it.
outname="output.tar.gz"
tar -czf ${outname} ${outdir_local}

# Ship the output tarball.
python copy_output.py -i ${outname} -e "tar.gz" -o ${outdir} -n ${proc_number}

# Cleanup. Not strictly necessary.
rm -r ${outdir_local}
rm ${outname}
rm *.py
rm -rf $gitdir
if [ -d fastjet ]; then rm -r fastjet; fi
if [ -d delphes ]; then rm -r delphes; fi
