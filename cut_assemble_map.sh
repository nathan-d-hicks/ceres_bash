#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=6-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=12  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=200G
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec
#SBATCH --job-name="megahit"

module load cutadapt
module load megahit
module load bowtie2
module load samtools

#Input is Read 1 name which has to have _R1_ in it like comes of casava

read1=$1
read2=`echo $read1 | sed 's/_R1_/_R2_/g'`
trim1=`echo $read1 | sed 's/.fastq.gz/.trimmed.fastq.gz/g'`
trim2=`echo $read2 | sed 's/.fastq.gz/.trimmed.fastq.gz/g'`
baseName=`echo $read1 | cut -d_ -f1`

cutadapt -j 24 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --max-n 3 -m 31 -o $trim1 -p $trim2 $read1 $read2 > ${baseName}.cutadapt.txt

myDir=`echo $PWD`
cd $TMPDIR

ln -s ${myDir}/${trim1} .
ln -s ${myDir}/${trim2} .

megahit -t 24 --k-min 31 --k-max 131 --k-step 10 -1 $trim1 -2 $trim2 -o megahit_out

cp megahit_out/log ${myDir}/${baseName}.log
cp megahit_out/options.json ${myDir}/${baseName}.options.json
cp megahit_out/final.contigs.fa ${myDir}/${baseName}.contigs.fa

cd $myDir

bowtie2-build ${baseName}.contigs.fa ${baseName}.contigs

bowtie2 -p 24 -x ${baseName}.contigs -1 $trim1 -2 $trim2 -S ${baseName}.sam 2> ${baseName}.bowtie.txt
samtools view -bS -F 4 ${baseName}.sam > ${baseName}.mapped.bam
samtools sort -o ${baseName}.S.mapped.bam ${baseName}.mapped.bam
samtools index ${baseName}.S.mapped.bam

rm ${baseName}.mapped.bam
rm ${baseName}.sam

echo "$baseName completed successfully"
