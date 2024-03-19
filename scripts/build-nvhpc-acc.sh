#!/bin/bash

module load nvhpc
# module load gcc/11.2.0-gcc-11.2.0
module load gcc/.12.3.0-gcc-11.2.0-nvptx
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-12.3.0-ab6j4u/lib64:$LD_LIBRARY_PATH
spack load netcdf-cxx4@4.3.1

# export TARGET="host"
export TARGET="multicore"
# export TARGET="gpu"

export CC=nvc
export CXX=nvc++
export CXXFLAGS="-std=c++17 -Minfo=accel -acc=$TARGET -gpu=managed"

export PREFIX="build-nvhpc-acc"

rm -rf $PREFIX

#configure muphys-cpp
cmake -DMU_IMPL=openacc -DMU_ENABLE_TESTS=OFF -B $PREFIX -S .

#build
cmake --build $PREFIX --parallel
