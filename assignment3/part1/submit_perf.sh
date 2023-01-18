#!/bin/sh
#BSUB -J proftest
#BSUB -q hpcintrogpu
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 10
#BSUB -R "rusage[mem=2048]"

export TMPDIR=$__LSF_JOB_TMPDIR__
module load cuda/11.8

export MFLOPS_MAX_IT=1

nv-nsight-cu-cli -o profile_$LSB_JOBID --set full ./matmult_c.nvc++ blk_offload 2048 2048 2048
