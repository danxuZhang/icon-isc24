#!/bin/bash

module purge
module load gcc/.12.3.0-gcc-11.2.0-nvptx
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-12.3.0-ab6j4u/lib64:$LD_LIBRARY_PATH
# module load gcc/.13.2.0-gcc-11.2.0-nvptx

spack load netcdf-cxx4@4.3.1

export CC=gcc
export CXX=g++
export CXXFLAGS="-std=c++17 -fopenacc -lm -foffload=nvptx-none -foffload-options=-lm"

export PREFIX="build-gcc-acc"

rm -rf $PREFIX

#configure muphys-cpp
cmake -DMU_IMPL=openacc -DMU_ENABLE_TESTS=OFF -B $PREFIX -S .

#build
cmake --build $PREFIX --parallel
