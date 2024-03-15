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

module load nvhpc
module load gcc/11.2.0-gcc-11.2.0

spack load netcdf-cxx4@4.3.1
spack load cdo@2.2.2

export INPUT="20k.nc"

export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64:$LD_LIBRARY_PATH
export ICON_DIR=/home/b/b382805/dx/icon-isc24
export SEQ_BUILD=$ICON_DIR/build-nvhpc-seq
export ACC_BUILD=$ICON_DIR/build-nvhpc-acc

cd $ICON_DIR

rm -f output.nc
rm -rf $SEQ_BUILD $ACC_BUILD

export CC=nvc
export CXX=nvc++
export CXXFLAGS=" -std=c++17 -Minfo=accel -fast -O3"

echo "Building Sequential: "
cmake -DMU_IMPL=seq -DMU_ENABLE_TESTS=OFF -B $SEQ_BUILD -S .
cmake --build $SEQ_BUILD --parallel 2>&1 | tee seq_build.log

echo "Building OpenACC: "
# export ACC_TARGET="multicore"
export ACC_TARGET="gpu"
export CXXFLAGS=" -acc=$ACCTARGET -gpu=managed $CXXFLAGS"
cmake -DMU_IMPL=openacc -DMU_ENABLE_TESTS=OFF -B $ACC_BUILD -S .
cmake --build $ACC_BUILD --parallel 2>&1 | tee acc_build.log

echo "Running Sequential:"
$SEQ_BUILD/bin/graupel ./tasks/$INPUT

echo "Evaulating Sequential Result: "
cdo infon -sub output.nc reference_results/sequential_double_$INPUT

echo "Running OpenACC: "
$ACC_BUILD/bin/graupel ./tasks/$INPUT
echo "Evaulating OpenACC Result: "
cdo infon -sub output.nc reference_results/sequential_double_$INPUT
