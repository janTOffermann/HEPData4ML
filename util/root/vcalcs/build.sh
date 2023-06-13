#!/bin/bash

# Build our custom ROOT library.
mkdir -p vectorcalcs/build
pushd vectorcalcs/build
  cmake ../../
  cmake --build .
  rm -r CMakeFiles
  rm **/*.cmake
popd
