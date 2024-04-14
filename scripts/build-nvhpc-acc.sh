#!/bin/bash

module load nvhpc/23.9-gcc-11.2.0

spack load netcdf-cxx4@4.3.1

# export TARGET="host"
# export TARGET="multicore"
export TARGET="gpu"

export CC=nvc
export CXX=nvc++
export CXXFLAGS="-std=c++17 -Minfo=accel -acc=$TARGET -gpu=managed -gpu=sm_80 -fast"

export PREFIX="build-nvhpc-acc"

rm -rf $PREFIX

#configure muphys-cpp
cmake -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_CXX_FLAGS="$CXXFLAGS" -DMU_IMPL=openacc -DMU_ENABLE_TESTS=OFF -B $PREFIX -S .

#build
cmake --build $PREFIX --parallel
