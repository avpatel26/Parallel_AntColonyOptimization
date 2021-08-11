#!/bin/bash
#SBATCH --partition=snowy
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=12:00:00

module load gcc/8.3.0
module load openmpi/3.1.4

srun  ./parallel.out graph.txt random.txt 32 15000 1
