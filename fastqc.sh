#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=10G
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec
#SBATCH --job-name="fastqc"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load fastqc

#make array of fastqs. Important if using wildcard to put quotes around it like "*gz"

readstring=`echo $1`

myOut="FASTQC_OUT"
mkdir $myOut

fastqc --noextract -t 16 -o $myOut $readstring
