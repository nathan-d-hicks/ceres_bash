#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=13-23:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=20G
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec
#SBATCH --job-name="bowtie-filter-sort"

module load bowtie2/2.3.4 
module load samtools

ref=$1
out=$2
bowtie2-build $1 $2

