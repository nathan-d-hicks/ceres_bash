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

#match string must be in single quotes if contains a wildcard
match_string=$1
#basename of the outfiles
outBase=$2

for read1 in $match_string
do
	read2=`echo $read1 | sed 's/_R1_/_R2_/g'`
	cat $read1 >> ${outBase}_1.fastq
	cat $read2 >> ${outBase}_2.fastq
done
