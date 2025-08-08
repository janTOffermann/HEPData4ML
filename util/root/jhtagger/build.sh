#!/bin/bash

# Add access to our CMake tools.
#export CMAKE_PREFIX_PATH={$CMAKE_PREFIX_PATH}:/local/home/jano/lgn/HEPData4ML/studies/W_subjet/cmake_util/cmaketools

# fj_libpath=$(realpath ../fastjet/fastjet-install/lib)
# fj_incpath=$(realpath ../fastjet/fastjet-install/include)
# CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:${fj_libpath}
# export FASTJET_INCLUDE_DIR=${fj_incpath}

rm -r jhtagger/build

# Now build our custom ROOT library.
mkdir -p jhtagger/build
pushd jhtagger/build
  cmake ../../
  cmake --build .
  rm -r CMakeFiles
  rm **/*.cmake
popd
