#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=6:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=2G
#SBATCH --partition=short
#SBATCH --job-name="compress"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


read1=$1
read2=`echo $read1 | sed 's/_R1_/_R2_/g'`
gzip -c $read1 > ${read1}.gz
gzip -c $read2 > ${read2}.gz



#stub=`echo $readName | cut -f1 -d"_"`
#mv  $readName.gz ${stub}.${readNum}.fastq.gz
