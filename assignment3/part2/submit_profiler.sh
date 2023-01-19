#!/bin/bash
#BSUB -J pois_prof
#BSUB -o pois_prof_%J.out
#BSUB -e pois_prof_%J.err
#BSUB -q hpcintrogpu
#BSUB -n 16
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 0:05
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2048]"

# The lock thing when profiling
export TMPDIR=$__LSF_JOB_TMPDIR__

module load nvhpc/22.11-nompi
module load cuda/11.8
module load gcc/11.3.0-binutils-2.38
make clean
make

command="poisson 500 1 0 0 -p gpu_mcp"

export MFLOPS_MAX_IT=1
#nsys profile -o profile_sys_$LSB_JOBID \
#    $command
nv-nsight-cu-cli -o profile_$LSB_JOBID \
    --section SpeedOfLight \
    --section Occupancy \
    --section MemoryWorkloadAnalysis \
    --section MemoryWorkloadAnalysis_Chart \
    --section ComputeWorkloadAnalysis \
    $command
