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
module load gcc/.12.3.0-gcc-11.2.0-nvptx
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-12.3.0-ab6j4u/lib64:$LD_LIBRARY_PATH

spack load netcdf-cxx4@4.3.1
spack load cdo@2.2.2

# export INPUT="20k.nc"
export INPUT="1500k.nc"

# export ACC_TARGET="multicore"
export ACC_TARGET="gpu"

export ICON_DIR=$SLURM_SUBMIT_DIR
export SEQ_BUILD=$ICON_DIR/build-nvhpc-seq
export ACC_BUILD=$ICON_DIR/build-nvhpc-acc

cd $ICON_DIR

rm -f output.nc seq_output.nc acc_output.nc
rm -rf $SEQ_BUILD $ACC_BUILD

export CC=nvc
export CXX=nvc++
export SEQ_FLAGS=" -O2 -std=c++17 -Wall -Wextra "
export ACC_FLAGS=" -acc=$ACC_TARGET -gpu=managed -gpu=sm_80 -gpu=fastmath -Minfo=accel -fast $SEQ_FLAGS"

echo "Building Sequential: "
cmake -DMU_IMPL=seq -DMU_ENABLE_TESTS=OFF \
    -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_CXX_FLAGS="$SEQ_FLAGS" \
    -DMU_ENABLE_SINGLE=OFF \
    -B $SEQ_BUILD -S .
cmake --build $SEQ_BUILD --parallel 2>&1 | tee seq_build.log

echo "Building OpenACC: "
cmake -DMU_IMPL=openacc -DMU_ENABLE_TESTS=OFF \
    -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_CXX_FLAGS="$ACC_FLAGS" \
    -DMU_ENABLE_SINGLE=OFF \
    -B $ACC_BUILD -S .
cmake --build $ACC_BUILD --parallel 2>&1 | tee acc_build.log

echo "Running Sequential:"
$SEQ_BUILD/bin/graupel ./tasks/$INPUT

mv output.nc seq_output.nc

echo "Running OpenACC: "
$ACC_BUILD/bin/graupel ./tasks/$INPUT
echo "Evaulating OpenACC Result: "
mv output.nc acc_output.nc
cdo infon -sub acc_output.nc seq_output.nc
# cdo infon -sub output.nc reference_results/sequential_double_$INPUT
