# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).

ROOT_STANDARD_LIBRARY_PACKAGE(NtupleProducer
  HEADERS
    ntuple/NtupleProducer.h
  SOURCES
    src/NtupleProducer.cxx
  LINKDEF
    NtupleProducer_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::RIO ROOT::MathCore ROOT::TreePlayer
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

target_include_directories(NtupleProducer PRIVATE
  ${HEPMC3_INCLUDE_DIR}
  "${CMAKE_SOURCE_DIR}/ntuple/inc"

)
target_link_libraries(NtupleProducer PRIVATE
  ${HEPMC3_LIBRARIES}
  ${HEPMC3_ROOT_LIBRARIES}
  ROOT::Core
  ROOT::RIO
  ROOT::MathCore
  ROOT::TreePlayer
)
