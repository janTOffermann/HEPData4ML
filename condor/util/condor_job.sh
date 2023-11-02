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
# $9a Event index offset.
# $10 Output directory (for the condor job).
# $11 Process number (for naming the output).
# $12 OpenBLAS max thread count (for multithreading).
# $13 Git option. Determines if we do a git clone here, or if the code has been shipped in as a tarball.
# $14 Git branch.

nevents_per_bin=$1
pt_bins=$2
separate_truth_flag=$3
separate_truth_number=$4
do_h5=$5
rng_seed=$6
do_split=$7
pythia_config=$8
event_idx_offset=$9
outdir=${10}
proc_number=${11}
openblas_max_thread=${12}
git_option=${13}
git_branch=${14}

# Set up the code. This may involve shipping in a payload, or running `git clone` here.
gitdir=HEPData4ML # TODO: payload curently set to use this name, is this OK or too much hardcoding?
if [ "${git_option}" == "0" ]; then
  # run git clone here
  git clone -b ${git_branch} git@github.com:janTOffermann/HEPData4ML.git ${gitdir}
else
  # assume the payload has been shipped in, as payload.tar.bz2
  payload=payload.tar.bz2
  tar -xjf ${payload}
  rm ${payload}
fi

# Run the setup script.
source ${gitdir}/setup/cvmfs/setup.sh

# Set the number of threads (for OpenBLAS). Might be necessary in order to deal with memory limits.
export OPENBLAS_NUM_THREADS=${openblas_max_thread}
export GOTO_NUM_THREADS=${openblas_max_thread}
export OMP_NUM_THREADS=${openblas_max_thread}

# Move the config.py file into the config directory. It has been shipped as an input file separate of the payload.
mv config.py ${gitdir}/config/config.py

outdir_local="output_${proc_number}"
truth_file_short=events.h5
delphes_file_short=events_delphes.h5
truth_file="${outdir_local}/${truth_file_short}"
delphes_file="${outdir_local}/${delphes_file_short}"

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

# Compress the output and extract it.
outname="output.tar.bz2"
tar -cjf ${outname} ${outdir_local}

# Ship the output tarball.
python copy_output.py -i ${outname} -e "tar.bz2" -o ${outdir} -n ${proc_number}

# Cleanup. Not strictly necessary.
rm -r ${outdir_local}
rm ${outname}
rm *.py
rm -rf $gitdir
if [ -d fastjet ]; then rm -r fastjet; fi
if [ -d delphes ]; then rm -r delphes; fi
