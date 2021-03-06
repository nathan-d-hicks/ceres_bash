#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=4  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=2G
#SBATCH --partition=short
#SBATCH --job-name="generic"

$1

