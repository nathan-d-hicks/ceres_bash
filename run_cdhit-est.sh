#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=24  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=10G
#SBATCH --job-name="cd-hit"
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load cd-hit

echo "Running:"
echo cd-hit-est $1 -i $2 -o $3
cd-hit-est $1 -T 24 -M 0 -i $2 -o $3


