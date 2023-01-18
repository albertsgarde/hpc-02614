#!/bin/bash
#BSUB -J proftest
#BSUB -q hpcintrogpu
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 10
#BSUB -R "rusage[mem=2048]"

bash analysis_2.sh