#!/bin/bash
#BSUB -J poisson
#BSUB -o poisson_%J.out
#BSUB -e poisson_%J.err
#BSUB -q hpcintrogpush
#BSUB -n 12
#BSUB -gpu "num=2:mode=exclusive_process"
#BSUB -W 0:05
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=512]"


module load nvhpc/22.11-nompi
module load cuda/11.8
module load gcc/11.3.0-binutils-2.38

NAME="dualGPU_test1.txt"

printf "Sequential\n" >> $NAME 
./poisson 200 200 0.1  -w -t -p seq >> $NAME
printf "With parallelization \n" >> $NAME
./poisson 200 200 0.1 -w -t -p par >> $NAME
printf "With Offloading using map\n"  >> $NAME
./poisson 200 200 0.1 -w -t -p gpu_map >> $NAME
printf "With Offloading and memcpy\n" >> $NAME
./poisson 200 200 0.1 -w -t -p gpu_mcp >> $NAME
printf "With Asynchrounous Offloading and dual GPUs  \n" >> $NAME
./poisson 200 200 0.1 -w2 -t -p gpu_async >> $NAME








