#!/bin/bash

eval $(/work/k20200/k202174/spack/bin/spack load --sh gcc@12.3.0)

spack load netcdf-cxx4@4.3.1

export CC=gcc
export CXX=g++
export CXXFLAGS="-std=c++17 -fopenacc -foffload=nvptx-none"

export PREFIX="build-gcc-acc"

rm -rf $PREFIX

#configure muphys-cpp
cmake -DMU_IMPL=openacc -DMU_ENABLE_TESTS=OFF -B $PREFIX -S .

#build
cmake --build $PREFIX --parallel
