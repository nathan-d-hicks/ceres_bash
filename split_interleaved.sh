#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=5G
#SBATCH --job-name="split_interleave"
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load bbtools


input=$1
R1out=`echo $input | sed 's/_filtered.fastq/_filtered.R1.fastq/g'`
R2out=`echo $input | sed 's/_filtered.fastq/_filtered.R2.fastq/g'`

time reformat.sh int=t addslash=t in=$input out1=$R1out out2=$R2out

