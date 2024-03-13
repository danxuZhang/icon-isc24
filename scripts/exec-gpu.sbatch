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

export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64:$LD_LIBRARY_PATH

export PREFIX=./build-nvhpc-acc

export INPUT="20k.nc"

cd /home/b/b382805/dx/icon-isc24

# echo "Profile run:"
module load nvhpc
nsys profile --trace=cuda,openacc --stats=true $PREFIX/bin/graupel ./tasks/20k.nc

$PREFIX/bin/graupel ./tasks/20k.nc
