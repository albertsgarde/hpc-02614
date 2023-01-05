#!/bin/bash
#BSUB -J mm_batch
#BSUB -o mm_batch_%J.out
#BSUB -q hpc
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 1:00
#BSUB -R "span[hosts=1] affinity[socket(1)]"

BLOCK_SIZES="1 4 16 64 256 1024"
bash analysis.sh ./matmult_c.gcc ./data_small "lib mnk mkn nmk nkm kmn knm" "10 50 100 250 375 500 750 1000 1250 1500" "$BLOCK_SIZES" ./cache_experiments_small
