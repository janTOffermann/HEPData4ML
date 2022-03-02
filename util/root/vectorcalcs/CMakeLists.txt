# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).
# Note that (as of October 2021), these ROOT CMake macros are not available if you've installed ROOT via conda.

ROOT_STANDARD_LIBRARY_PACKAGE(VectorCalcs
  HEADERS
    vectorcalcs/VectorCalcs.h
  SOURCES
    src/VectorCalcs.cxx
  LINKDEF
    VectorCalcs_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::Tree ROOT::RIO ROOT::MathCore
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM

)
