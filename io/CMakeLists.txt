# find netcdf
find_package(NetCDF COMPONENTS C CXX REQUIRED)
include_directories(${NetCDF_CXX_INCLUDE_DIR})

add_library(muphys_io SHARED "io.cpp")
set_target_properties(muphys_io PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(muphys_io PUBLIC NetCDF::NetCDF_CXX NetCDF::NetCDF_C)
