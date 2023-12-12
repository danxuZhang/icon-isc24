#!/bin/bash
  
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

ulimit -s unlimited
ulimit -c 0

# compiler = gnu/intel/nvidia
COMPILER='gnu' 

# compiler flags
FLAGS="" 

# clean or dirty build
FORCE_CLEAN=true

. scripts/levante-setup.sh $COMPILER gpu

# TODO: define scripts