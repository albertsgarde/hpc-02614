#!/bin/bash
#BSUB -J poisson
#BSUB -o poisson_%J.out
#BSUB -e poisson_%J.err
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=128]"
#BSUB -W 1:00
#BSUB -R "span[hosts=1] affinity[socket(1)]"

module load gcc/12.2.0-binutils-2.39
make realclean
make
. ../../env/bin/activate
python3 experiments/parallel_running_time.py

