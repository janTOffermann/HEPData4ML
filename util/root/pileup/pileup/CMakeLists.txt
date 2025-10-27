# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).

ROOT_STANDARD_LIBRARY_PACKAGE(Pileup
  HEADERS
    pileup/Pileup.h
  SOURCES
    src/Pileup.cxx
  LINKDEF
    Pileup_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::RIO ROOT::MathCore ROOT::Hist
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

target_include_directories(Pileup PRIVATE
  "${HEPMC3_INCLUDE_DIR}"
    "${CMAKE_SOURCE_DIR}/pileup/inc"

)
target_link_libraries(Pileup PRIVATE
  "${HEPMC3_LIBRARIES}"
  "${HEPMC3_ROOT_LIBRARIES}"
  ROOT::Core
  ROOT::RIO
  ROOT::MathCore
  ROOT::Hist
)
