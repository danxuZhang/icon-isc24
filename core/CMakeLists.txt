add_library(muphys_core SHARED "common/utils.cpp" "common/graupel.hpp")
target_include_directories(muphys_core PUBLIC common properties transitions)
set_target_properties(muphys_core PROPERTIES LINKER_LANGUAGE CXX)