#!/bin/bash

#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.log
#SBATCH --time=00:10:00

ulimit -s unlimited
ulimit -c 0

# compiler = gnu/intel/nvidia
COMPILER='gnu' 

# compiler flags
FLAGS="" 

# clean or dirty build
FORCE_CLEAN=true

# input file
FILE=$(pwd)/tasks/input.nc

. scripts/levante-setup.sh $COMPILER cpu

# build the code
. scripts/build.sh "build" "-DMU_IMPL=seq -DCMAKE_CXX_FLAGS=$FLAGS" $FORCE_CLEAN

# run the executable
./build/bin/graupel $FILE