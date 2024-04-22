#!/bin/bash

#SBATCH --account=ka1273_gpu
#SBATCH --job-name=icon-single-percision
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

ulimit -s unlimited
ulimit -c 0

module purge
module load nvhpc/23.9-gcc-11.2.0
module load gcc/.12.3.0-gcc-11.2.0-nvptx
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-12.3.0-ab6j4u/lib64:$LD_LIBRARY_PATH

spack load netcdf-cxx4@4.3.1
spack load cdo@2.2.2

# export INPUT="20k.nc"
export INPUT="1500k.nc"

# export ACC_TARGET="multicore"
export ACC_TARGET="gpu"

export ICON_DIR=$SLURM_SUBMIT_DIR
export BASE_BUILD=$ICON_DIR/build-base
export SEQ_BUILD=$ICON_DIR/build-gcc-seq
export ACC_BUILD=$ICON_DIR/build-nvhpc-acc

cd $ICON_DIR

rm -f output.nc seq_output.nc acc_output.nc
rm -rf $SEQ_BUILD $ACC_BUILD $BASE_BUILD

export CC=nvc
export CXX=nvc++
export SEQ_FLAGS=" -O3 -std=c++17 -Wall -Wextra "
export ACC_FLAGS=" -acc=$ACC_TARGET -gpu=managed -gpu=cc80 -gpu=fastmath -Minfo=accel -fast $SEQ_FLAGS"

echo "========================================================="
echo "Building OpenACC impl with gcc for correctness checking: "
cmake -DMU_IMPL=openacc -DMU_ENABLE_TESTS=ON \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_CXX_FLAGS="-O0" \
    -DMU_ENABLE_SINGLE=ON \
    -B $SEQ_BUILD -S .
cmake --build $SEQ_BUILD --parallel 2>&1 | tee seq_build.log
echo "========================================================="

echo ""
echo ""

echo "========================================================="
echo "Build Sequential impl with nvhpc for comparison baseline:"
cmake -DMU_IMPL=seq -DMU_ENABLE_TESTS=OFF \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_CXX_FLAGS="-O2" \
    -DMU_ENABLE_SINGLE=ON \
    -B $BASE_BUILD -S .
cmake --build $BASE_BUILD --parallel 2>&1 | tee base_build.log
echo "========================================================="

echo ""
echo ""

echo "========================================================="
echo "Building OpenACC impl: "
cmake -DMU_IMPL=openacc -DMU_ENABLE_TESTS=OFF \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_CXX_FLAGS="$ACC_FLAGS" \
    -DMU_ENABLE_SINGLE=ON \
    -B $ACC_BUILD -S .
cmake --build $ACC_BUILD --parallel 2>&1 | tee acc_build.log

echo "========================================================="
echo "Running OpenACC impl gcc build 20k:"
$SEQ_BUILD/bin/graupel ./tasks/20k.nc
echo "Comparing with reference result:"
cdo diffn output.nc ./reference_results/sequential_single_20k.nc
echo "Return status: $?"
mv output.nc acc_single_20k.nc
echo "========================================================="

echo ""
echo ""

echo "========================================================="
echo "Running Sequential impl (baseline) 20k:"
$BASE_BUILD/bin/graupel ./tasks/20k.nc
rm output.nc

echo "Running OpenACC impl: 20k"
$ACC_BUILD/bin/graupel ./tasks/20k.nc
echo "Compare OpenACC GPU output with Reference: "
cdo infon -sub output.nc ./reference_results/sequential_single_20k.nc
rm output.nc
echo "========================================================="

echo ""
echo ""

echo "========================================================="
echo "Running Sequential impl (baseline) 1500k:"
$BASE_BUILD/bin/graupel ./tasks/1500k.nc
mv output.nc base_output.nc

echo "Running OpenACC impl: 1500k "
$ACC_BUILD/bin/graupel ./tasks/1500k.nc
echo "Evaulating OpenACC Result: "
mv output.nc acc_output.nc
echo "Compare OpenACC GPU output with Baseline output: "
cdo infon -sub acc_output.nc base_output.nc
echo "========================================================="

# rm acc_output.nc seq_output.nc
# nsys profile --trace=cuda,openacc --stats=true --gpu-metrics-device all -o 20k-128-128-128-exp3 $ACC_BUILD/bin/graupel ./tasks/20k.nc
# rm output.nc
# nsys profile --trace=cuda,openacc --stats=true --gpu-metrics-device all -o 1500k-128-128-128-exp3 $ACC_BUILD/bin/graupel ./tasks/1500k.nc
# rm output.nc
