#!/bin/bash

module load nvhpc
module load gcc/11.2.0-gcc-11.2.0

spack load netcdf-cxx4@4.3.1

export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64:$LD_LIBRARY_PATH

export CC=nvc
export CXX=nvc++
export CXXFLAGS="-std=c++17 -Minfo=all"

export PREFIX="build-nvhpc-seq"

#configure muphys-cpp
cmake -DMU_IMPL=seq -B $PREFIX -S .

#build 
cmake --build $PREFIX --parallel
