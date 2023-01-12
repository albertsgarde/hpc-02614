#!/bin/bash
#BSUB -J poisson
#BSUB -o poisson_%J.out
#BSUB -e poisson_%J.err
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -R "rusage[mem=2048]"
#BSUB -W 0:30

module load gcc/12.2.0-binutils-2.39
make realclean
make
. ../../env/bin/activate
python3 experiments/parallel.py