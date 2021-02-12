#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=4  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=4G
#SBATCH --partition=short
#SBATCH --job-name="cutadapt"

module load bbtools

#Input is Read 1 name which has to have _R1_ in it like comes of casava. Need to specify trimmed reads

read1=$1
read2=`echo $read1 | sed 's/_R1_/_R2_/g'`
baseName=`echo $read1 | cut -d_ -f1`
out1=${baseName}.R1.fastq.gz
out2=${baseName}.R2.fastq.gz
outMerge=${baseName}.merge.fastq.gz

bbmerge.sh in1=$read1 in2=$read2 out=$outMerge outu1=$out1 outu2=$out2
