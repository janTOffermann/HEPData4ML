# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).

ROOT_STANDARD_LIBRARY_PACKAGE(JHTagger
  HEADERS
    jhtagger/JHTagger.h
  SOURCES
    src/JHTagger.cxx
  LINKDEF
    JHTagger_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::RIO ROOT::MathCore
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

# Now, let's try to also link the fastjet stuff -- this isn't using the above ROOT macro,
# so hopefully it still works (though it makes this all look a bit ugly).
# It requires that the FASTJET_INCLUDE_DIR and FASTJET_LIB variables are set.
target_include_directories(JHTagger PRIVATE "${FASTJET_INCLUDE_DIR}")
target_link_libraries(JHTagger PRIVATE "${FASTJET_LIB}")
