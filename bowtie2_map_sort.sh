#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=10G
#SBATCH --partition=short
#SBATCH --job-name="bowtie-filter-sort"

module load bowtie2/2.3.4 
module load samtools


read1=$1
read2=$2
outStub=$3
ref=$4

bowtie2 -p 20 -x $ref -1 $read1 -2 $read2 -S ${outStub}.sam 2> ${outStub}.bowtie.txt
samtools view -b -F 4 ${outStub}.sam > ${outStub}.mapped.bam
samtools sort -o ${outStub}.S.mapped.bam ${outStub}.mapped.bam
samtools index ${outStub}.S.mapped.bam

rm ${outStub}.sam
rm ${outStub}.mapped.bam
