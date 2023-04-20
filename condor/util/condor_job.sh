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
# $9 Output directory (for the condor job).
# $10 Process number (for naming the output).

# Unpack the tarball that contains our generation code.
tar -xjf payload.tar.bz2
mv payload/* ./
rm payload.tar.bz2
rm -r payload

# Run the setup script.
source setup/cvmfs/setup.sh

# Delete any VectorCalcs build that might've been included in the payload,
# it will have to be rebuilt locally based on how it's currently implemented.
# (the build is very quick)
vc_build=util/root/vectorcalcs/build
if [ -d $vc_build ]; then rm -r $vc_build; fi

# Move the config.py file into the config directory. It has been shipped as an input file separate of the payload.
mv config.py config/config.py

outdir_local="output_${10}"

# Run the generation.
python run.py -n $1 -p_s $2 -O ${outdir_local} -s $3 -ns $4 -h5 $5 -rng $6 -pb 1 --split $7 -cd 1 -pc $8 # will have to put in a bunch of user options.

# Compress the output and extract it.
outname="output.tar.bz2"
tar -cjf ${outname} ${outdir_local}
rm -r ${outdir_local}

# Ship the output.
python copy_output.py -i ${outname} -e "tar.bz2" -o ${9} -n ${10}

# Cleanup. Not strictly necessary.
rm ${outname}
rm *.py
rm -r util
rm -r config
rm -r setup
if [ -d fastjet ]; then rm -r fastjet; fi
if [ -d delphes ]; then rm -r delphes; fi