#!/bin/bash

module load nvhpc/22.11-nompi
module load cuda/11.8
module load gcc/11.3.0-binutils-2.38

# printf "Sequential\n"
# ./poisson 30 50 0.1 -t -p seq
# printf "With parallelization \n"
# ./poisson 30 200 0.1 -t -p par
# printf "With Offloading using map\n"
# ./poisson 30 200 0.1 -t -p gpu_map
# printf "With Offloading and memcpy\n"
# ./poisson 30 200 0.1 -t -p gpu_mcp
printf "With Asynchrounous Offloading  \n"
./poisson 10 10 0.1 -t -p gpu_async
./poisson 10 100 0.1 -t -p gpu_async


