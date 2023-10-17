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

# Run the generation.
python ${gitdir}/run.py -n $1 -p_s $2 -O ${outdir_local} -s $3 -ns $4 -h5 $5 -rng $6 -pb 1 --split $7 -cd 1 -pc $8 -df 1 --index_offset $9 # will have to put in a bunch of user options.

# Compress the output and extract it.
outname="output.tar.bz2"
tar -cjf ${outname} ${outdir_local}

# Ship the output tarball.
python copy_output.py -i ${outname} -e "tar.bz2" -o ${10} -n ${11}

# Additionally, ship the full event HDF5 file(s).
# This is technically redundant but might save us some time, since we don't need to unpack all the job output for them.
if [ "${7}" == "1" ]; then
  python copy_output.py -i ${outdir_local}/train.h5 -e "h5" -o ${10} -n ${11}
  python copy_output.py -i ${outdir_local}/test.h5  -e "h5" -o ${10} -n ${11}
  python copy_output.py -i ${outdir_local}/valid.h5 -e "h5" -o ${10} -n ${11}
else
  python copy_output.py -i ${outdir_local}/events.h5 -e "h5" -o ${10} -n ${11}
fi

# Cleanup. Not strictly necessary.
rm -r ${outdir_local}
rm ${outname}
rm *.py
rm -r $gitdir
if [ -d fastjet ]; then rm -r fastjet; fi
if [ -d delphes ]; then rm -r delphes; fi
