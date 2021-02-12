#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=15G
#SBATCH --job-name="prokka"
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load prokka

fasta=$1
out=`echo $fasta | sed 's/\.[^.]*$//'`"_prokka"

prokka --cpus 8 --prefix $out $fasta 


