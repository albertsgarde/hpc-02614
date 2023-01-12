#!/bin/bash
#BSUB -J poisson
#BSUB -o poisson_%J.out
#BSUB -e poisson_%J.err
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -R "rusage[mem=512]"
#BSUB -W 1:00
#BSUB -R "span[hosts=1]"


module load gcc/12.2.0-binutils-2.39
make realclean
make
echo "Activating python environment."
. ../../env/bin/activate
echo "Running experiment script."
python3 experiments/parallel_running_time.py

