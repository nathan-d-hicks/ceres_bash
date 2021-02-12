#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=4  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --job-name="cutadapt"

module load miniconda
source activate cutadaptenv


#Input is Read 1 name which has to have _R1_ in it like comes of casava

read1=$1
read2=`echo $read1 | sed 's/_R1_/_R2_/g'`
out1=`echo $read1 | sed 's/.fastq.gz/.trimmed.fastq.gz/g'`
out2=`echo $read2 | sed 's/.fastq.gz/.trimmed.fastq.gz/g'`
outLog=`echo $read1 | cut -d_ -f1`

cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 5 -m 31 -o $out1 -p $out2 $read1 $read2 > ${outLog}.cutadapt.txt
