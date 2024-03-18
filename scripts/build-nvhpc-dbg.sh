#!/bin/bash

module load nvhpc
module load gcc/.12.3.0-gcc-11.2.0-nvptx
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-12.3.0-ab6j4u/lib64:$LD_LIBRARY_PATH

spack load netcdf-cxx4@4.3.1

export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64:$LD_LIBRARY_PATH

export CC=nvc
export CXX=nvc++
export CXXFLAGS=" -std=c++17 -g "

# export IMPLEMENTATION="seq"
export IMPLEMENTATION="openacc"

export PREFIX="build-nvhpc-dbg"

rm -rf $PREFIX
#configure muphys-cpp
cmake -DMU_IMPL=$IMPLEMENTATION -DCMAKE_BUILD_TYPE=Debug -DMU_ENABLE_TESTS=OFF -B $PREFIX -S .

#build
cmake --build $PREFIX --parallel
