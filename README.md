# muphys-cpp

## Description
C++ prototypes of [muphys](https://gitlab.dkrz.de/icon-libraries/muphys) project using heterogeneous libraries and C extensions.

## Documentation
Online documentation is generated automatically using `doxygen` .

## Installation

### Dependencies
* [NetCDF for CXX](https://github.com/Unidata/netcdf-cxx4)
  * for Levante: `spack load netcdf-cxx4@4.3.1`

Other dependency like [googletest](https://github.com/google/googletest) is built in-tree from github archives. 

### Available compile options 
* _Implementation_ - The sequential implementation is selected by default. The user can choose of the following options:
  * MU_IMPL=seq - C++ serial implementation
 * _Precision_ (default is `double`)
  * MU_ENABLE_SINGLE - to switch to `float` 
* _Unit-test_ - compile tests together with the main executable (default is `true`)
  * MU_ENABLE_TESTS

### Compile the project (with default flags and Seq frontend)

`cmake -DMU_IMPL=seq -B build -S .`

`cmake --build build`

## Usage

`./build/bin/graupel input_file.nc`

### Automated tests

- Run tests manually:
`cd build && ctest` 

## License

muphys-cpp is available under a BSD 3-clause license. See [LICENSES/](./LICENSES) for license information and [AUTHORS.TXT](./AUTHORS.TXT) for a list of authors.
