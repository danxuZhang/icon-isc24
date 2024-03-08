#!/bin/bash
  
#SBATCH --account=ka1273_gpu
#SBATCH --job-name=icon
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

ulimit -s unlimited
ulimit -c 0

export PREFIX=./build-nvhpc-gpu

cd /home/b/b382805/dx/icon-isc24
$PREFIX/bin/graupel ./tasks/20k.nc
