#!/bin/bash

buildDir="display/build"

# if [ -d $buildDir ]; then
#   rm -r $buildDir
# fi

# Now build our custom ROOT library.
mkdir -p $buildDir
pushd $buildDir
  cmake ../../
  cmake --build .
  rm -r CMakeFiles
  rm **/*.cmake
popd
