#!/bin/bash
#BSUB -J jacobi_par
#BSUB -o jacobi_par_%J.out
#BSUB -q hpcintro
#BSUB -n 16
#BSUB -R "rusage[mem=2048]"
#BSUB -W 0:10

rm parallel_j_results.dat

NS="1 2 4 8 16"
for n in $NS
do
  OMP_NUM_THREADS=$n ./poisson_j 100 10000 0 0.0 >> parallel_j_results.dat
done
