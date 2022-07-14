#!/bin/bash
# Set up the environment (make sure conda environment is active).
# conda activate lgn_data

# Now build our custom ROOT library.
mkdir -p vectorcalcs/build
pushd vectorcalcs/build
  cmake ../../
  cmake --build .
  rm -r CMakeFiles
  rm **/*.cmake
popd
