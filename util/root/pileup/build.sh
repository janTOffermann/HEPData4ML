#!/bin/bash
buildDir="pileup/build"

rm -r $buildDir

# Now build our custom ROOT library.
mkdir -p $buildDir
pushd $buildDir
  cmake ../../
  cmake --build .
  rm -r CMakeFiles
  rm **/*.cmake
popd
