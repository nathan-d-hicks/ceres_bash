#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=100G
#SBATCH --job-name="sortmerna"
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load sortmerna

r5p8S='/project/epicon/bin/rRNA_databases/rfam-5.8s-database-id98.fasta'
r5S='/project/epicon/bin/rRNA_databases/rfam-5s-database-id98.fasta'
r16S='/project/epicon/bin/rRNA_databases/silva-bac-16s-id90.fasta'
rarc16S='/project/epicon/bin/rRNA_databases/silva-arc-16s-id95.fasta'
r18S='/project/epicon/bin/rRNA_databases/silva-euk-18s-id95.fasta'
r23S='/project/epicon/bin/rRNA_databases/silva-bac-23s-id98.fasta'
rarc23S='/project/epicon/bin/rRNA_databases/silva-arc-23s-id98.fasta'
r28S='/project/epicon/bin/rRNA_databases/silva-euk-28s-id98.fasta'

read1=$1
read2=`echo $read1 | sed 's/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz/g'`
outstub=`echo $read1 | sed 's/_1.trimmed.fastq.gz//g'`

time sortmerna -ref $r5p8S -ref $r5S -ref $r16S -ref $rarc16S -ref $r18S -ref $r23S -ref $rarc23S -ref $r28S -reads $read1 -reads $read2 -fastx -other -paired_in -workdir ${outstub}_sortmerna -a 8

cp ${outstub}_sortmerna/out/aligned.fastq ${outstub}_rRNA.fastq
cp ${outstub}_sortmerna/out/other.fastq ${outstub}_filtered.fastq
