#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=6-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=12  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=180G
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
baseName=`echo $read1 | cut -d_ -f1`
echo $baseName
trim1=${baseName}.R1.trim.fastq.gz
trim2=${baseName}.R2.trim.fastq.gz
echo $trim1
echo $trim2

#trim off the adapters
cutadapt -j 24 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --max-n 3 -m 31 -o $trim1 -p $trim2 $read1 $read2 > ${baseName}.cutadapt.txt

#map to the sorghum genome
bowtie2 -p 24 -x /project/epicon/bin/sorghum_genome/Sbicolor_454_v3.0.1.hardmasked -1 $trim1 -2 $trim2 -S ${baseName}.sam 2> ${baseName}.sorghum.bowtie.txt

#remove reads where mapped or mate NOT  mapped (f 12) and where not primary alignemnt (F 256)
samtools view -bS -f 12 -F 256 ${baseName}.sam > ${baseName}.sorghum.bam
rm ${baseName}.sam
samtools sort -n -m 5G -@ 8 ${baseName}.sorghum.bam -o ${baseName}.sorghum.S.bam

#now write out the only paired ones into -1 and -2 so that this is where neither read mapped to the sorghum genome
samtools fastq -@ 8 ${baseName}.sorghum.S.bam -1 ${baseName}.R1.filter.fastq.gz -2 ${baseName}.R2.filter.fastq.gz -s /dev/null

rm ${baseName}.sorghum.bam
rm ${baseName}.sorghum.S.bam

myDir=`echo $PWD`
cd $TMPDIR

ln -s ${myDir}/${baseName}.R1.filter.fastq.gz .
ln -s ${myDir}/${baseName}.R2.filter.fastq.gz .

megahit -t 24 --k-min 31 --k-max 131 --k-step 10 -1 ${baseName}.R1.filter.fastq.gz -2 ${baseName}.R2.filter.fastq.gz -o ${baseName}_megahit

cp ${baseName}_megahit/log ${baseName}.log
cp ${baseName}_megahit/options.json ${baseName}.options.json
cp ${baseName}_megahit/final.contigs.fa ${baseName}.contigs.fa

cd $myDir

bowtie2-build ${baseName}.contigs.fa ${baseName}.contigs

bowtie2 -p 24 -x ${baseName}.contigs -1 ${baseName}.R1.filter.fastq.gz -2 ${baseName}.R2.filter.fastq.gz -S ${baseName}.sam 2> ${baseName}.self.bowtie.txt
samtools view -bS -F 4 ${baseName}.sam > ${baseName}.mapped.bam
samtools sort -o ${baseName}.S.mapped.bam ${baseName}.mapped.bam
samtools index ${baseName}.S.mapped.bam

rm ${baseName}.mapped.bam
rm ${baseName}.sam

echo "$baseName completed successfully"
