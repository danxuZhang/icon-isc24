include(FetchContent)

enable_testing()

# include all source files
include_directories(${CMAKE_SOURCE_DIR})
   
# create test executable
add_executable(muphys_core_test common.cc io.cc)

target_include_directories(muphys_core_test PRIVATE ${CMAKE_SOURCE_DIR})

# link against googletest (built locally from /extern)
target_link_libraries(muphys_core_test GTest::gtest_main muphys_core muphys_io muphys_implementation)

include(GoogleTest)
gtest_discover_tests(muphys_core_test)

add_test(
  NAME tiny
  COMMAND ${CMAKE_BINARY_DIR}/bin/graupel ${CMAKE_BINARY_DIR}/test/input/tiny.nc
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
FetchContent_Declare(tiny
  URL "https://swift.dkrz.de/v1/dkrz_ab6243f85fe24767bb1508712d1eb504/muphys_data/tiny.nc?temp_url_sig=edb55a3be02227d9606eec8ec06fc066cb105110&temp_url_expires=2026-09-03T09:09:43Z"
  URL_MD5 2996eeba471abf2a614257ea80939dc9
  DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/test/input
  DOWNLOAD_NAME "tiny.nc"
  DOWNLOAD_NO_EXTRACT True
)
FetchContent_MakeAvailable(tiny)

add_test(
  NAME mini
  COMMAND ${CMAKE_BINARY_DIR}/bin/graupel ${CMAKE_BINARY_DIR}/test/input/mini.nc
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
FetchContent_Declare(mini
  URL "https://swift.dkrz.de/v1/dkrz_ab6243f85fe24767bb1508712d1eb504/muphys_data/atm_R2B03.nc?temp_url_sig=03a3ddf3d07e1b40c609eb0a24b3a9f78cda689d&temp_url_expires=2026-08-05T18:40:02Z"
  URL_MD5 3b4db7ed1ca361bba8111701eb616cde
  DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/test/input
  DOWNLOAD_NAME "mini.nc"
  DOWNLOAD_NO_EXTRACT True
)
FetchContent_MakeAvailable(mini)

add_test(
  NAME dbg
  COMMAND ${CMAKE_BINARY_DIR}/bin/graupel ${CMAKE_BINARY_DIR}/test/input/dbg.nc
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
FetchContent_Declare(dbg
  URL "https://swift.dkrz.de/v1/dkrz_ab6243f85fe24767bb1508712d1eb504/muphys_data/dbg.nc?temp_url_sig=2ffaaed83d57b3501e7fd6a343c13a3e64af931a&temp_url_expires=2026-08-03T14:45:09Z"
  URL_MD5 39422f4a9f35b70e00107bedb6644c8f
  DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/test/input
  DOWNLOAD_NAME "dbg.nc"
  DOWNLOAD_NO_EXTRACT True
)
FetchContent_MakeAvailable(dbg)
