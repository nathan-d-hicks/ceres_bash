#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=2-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=12  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=10G
#SBATCH --job-name="kallisto"
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# Run as run_kallisto.sh <R1> <R2> <output_folder>

module load miniconda
source activate ndh_env

#kallisto index -i Prokka_rhizo_coassembled $1

kallisto quant -i 55fa_genes.k.index --single-overhang --rf-stranded -t 24 -o $3 $1 $2
