#!/bin/bash
#BSUB -J poisson
#BSUB -o poisson_%J.out
#BSUB -e poisson_%J.err
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=128]"
#BSUB -W 0:30
#BSUB -R "span[hosts=1] affinity[socket(1)]"

make realclean
make
. ../../env/bin/activate
python3 experiments/convergence.py

