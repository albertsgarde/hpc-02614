#!/bin/bash
#BSUB -J cache_analyses_batch
#BSUB -o cache_analyses_batch_%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 60
#BSUB -R "span[hosts=1] affinity[socket(1)]"

bash cache_analysis.sh