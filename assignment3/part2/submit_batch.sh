#!/bin/bash
#BSUB -J poisson
#BSUB -o poisson_%J.out
#BSUB -e poisson_%J.err
#BSUB -q hpcintrogpu
#BSUB -n 16
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 0:05
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=512]"

# The lock thing when profiling
export TMPDIR=$__LSF_JOB_TMPDIR__

module load nvhpc/22.11-nompi
module load cuda/11.8
module load gcc/11.3.0-binutils-2.38
make clean
make
echo "Activating python environment."
. ../../env/bin/activate
echo "Running experiment script."
python3 experiments/frobenius_cost.py
