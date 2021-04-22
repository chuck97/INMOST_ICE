#!/bin/bash
#SBATCH --job-name=CG2_FCT
#SBATCH --ntasks=32
#SBATCH --time=12:00:00
#SBATCH --partition=mix

mpirun ./run ./config.json > ./output.txt
