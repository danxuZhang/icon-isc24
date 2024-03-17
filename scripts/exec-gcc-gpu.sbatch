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

module load gcc/.12.3.0-gcc-11.2.0-nvptx
export PREFIX=./build-gcc-acc

export INPUT="20k.nc"

# echo "Profile run:"
module load nvhpc
nsys profile --trace=cuda,openacc --stats=true $PREFIX/bin/graupel ./tasks/20k.nc

echo "Running graupel"
cd /home/b/b382805/dx/icon-isc24
$PREFIX/bin/graupel ./tasks/$INPUT
