#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=13-23:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=30   # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=800G
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec
#SBATCH --job-name="megahit"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#usage: 

#sbatch ~/bin/megahit.sh '*.R1.trimmed.fastq.gz’

module load megahit/1.2.9

cd $TMPDIR

#input is a list of fastq including as a wildcard
#must have the .R1. .R2 style as I cahnge them to
#link them all into the tmp directory on the node

for read1 in $1
do
	read2=`echo $read1 | sed 's/.R1./.R2./g'`
	ln -s /project/epicon/nhicks/${read1} .
	ln -s /project/epicon/nhicks/${read2} .
done

# this is a one liner I don't totally understand but will join together stuff with a single character that you specify
function join_by { local IFS="$1"; shift; echo "$*"; }
read1list=`join_by , $1`
read2list=`echo $read1list | sed 's/.R1./.R2./g'`

megahit -t 30 --k-min 31 --k-max 201 --k-step 10 -1 $read1list -2 $read2list -o megahit_out

baseName=`echo $read1 | cut -d. -f1`

cp megahit_out/log ${myDir}/${baseName}.log
cp megahit_out/options.json ${myDir}/${baseName}.options.json
cp megahit_out/final.contigs.fa ${myDir}/${baseName}.contigs.fa


