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

# Clone the git repo -- this is cleaner than packaging up the code with a tarball.
git clone -b devel git@github.com:janTOffermann/HEPData4ML.git # get the development branch, maybe make this configurable?
gitdir=HEPData4ML

# Run the setup script.
source ${gitdir}/setup/cvmfs/setup.sh

# Set the number of threads (for OpenBLAS). Might be necessary in order to deal with memory limits.
export OPENBLAS_NUM_THREADS=$12
export GOTO_NUM_THREADS=$12
export OMP_NUM_THREADS=$12

# Move the config.py file into the config directory. It has been shipped as an input file separate of the payload.
mv config.py ${gitdir}/config/config.py

outdir_local="output_${11}"

truth_file="${outdir_local}/events.h5"
delphes_file="${outdir_local}/event_delphes.h5"

# Run the generation. First with Delphes, then again without.
# (for the second run, use the HepMC files that have already been generated)
python ${gitdir}/run.py \
  -n $1 \
  -p_s $2 \
  -O ${outdir_local} \
  -s $3 \
  -ns $4 \
  -h5 $5 \
  -rng $6 \
  -pb 1 \
  --delete_stats 0 \
  --split 0 \
  -del_delphes 1 \
  -pc $8 \
  -df 1 \
  --index_offset $9 \
  --delphes 1

mv $truth_file $delphes_file # first output is from Delphes -> rename it.

# Run a 2nd time, now without Delphes output. Produces truth-level "events.h5"
python ${gitdir}/run.py \
  -n $1 \
  -p_s $2 \
  -O ${outdir_local} \
  -s $3 \
  -ns $4 \
  -h5 $5 \
  -rng $6 \
  -pb 1 \
  --delete_stats 0 \
  --split 0 \
  -del_delphes 1 \
  -pc $8 \
  -df 1 \
  --index_offset $9 \
  --delphes 0

# TODO: Optionally do something here with events.h5 and events_delphes.h5 -- slim down to only shared events
#       that passed both truth-level and Delphes jet selections?

# Optionally split files into train,test and validation sets.
# TODO: Split is currently hard-coded, make this configurable?
if [ "${7}" == "1" ]; then
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
if [ "${7}" == "1" ]; then
  python copy_output.py -i ${outdir_local}/train.h5 -e "h5" -o ${10} -n ${11}
  python copy_output.py -i ${outdir_local}/test.h5  -e "h5" -o ${10} -n ${11}
  python copy_output.py -i ${outdir_local}/valid.h5 -e "h5" -o ${10} -n ${11}
  rm ${outdir_local}/train.h5 ${outdir_local}/test.h5 ${outdir_local}/valid.h5

  python copy_output.py -i ${outdir_local}/delphes_train.h5 -e "h5" -o ${10} -n ${11}
  python copy_output.py -i ${outdir_local}/delphes_test.h5  -e "h5" -o ${10} -n ${11}
  python copy_output.py -i ${outdir_local}/delphes_valid.h5 -e "h5" -o ${10} -n ${11}
  rm ${outdir_local}/delphes_train.h5 ${outdir_local}/delphes_test.h5 ${outdir_local}/delphes_valid.h5

else
  python copy_output.py -i $truth_file -e "h5" -o ${10} -n ${11}
  python copy_output.py -i $delphes_file -e "h5" -o ${10} -n ${11}

  rm $truth_file
  rm $delphes_file
fi

# Compress the output and extract it.
outname="output.tar.bz2"
tar -cjf ${outname} ${outdir_local}

# Ship the output tarball.
python copy_output.py -i ${outname} -e "tar.bz2" -o ${10} -n ${11}

# Cleanup. Not strictly necessary.
rm -r ${outdir_local}
rm ${outname}
rm *.py
rm -rf $gitdir
if [ -d fastjet ]; then rm -r fastjet; fi
if [ -d delphes ]; then rm -r delphes; fi
