# 16s Qiime Analysis with DADA2






# WGS Metagenome Analysis

## Getting data, manipulating fastq files

**submit_download_sra.sh**: Control file that submits download_sra.sh for every SRA number that is listed one per line in a plain text file.  So if you have 20 lines of SRA numbers, you get 20 slurm jobs downloading this.

`submit_download_sra.sh <sra_list_file>`

**download_sra.sh**: uses fasterq-dump to download SRA data as fastq

```
downlaod_sra.sh SRA#
```

**combine_reads.sh**: takes a wildcard matching read 1 data and combines all the matching fastQ into a single fastq file.  Set to take unzipped now, but just need to change the name if the files are .gz since you can concatenate .gz files.  Yields outstub_1.fastq and outstub_2.fastq

```combine_reads.sh '*_R1_*' <out_stub>```


## QC and Pre-processing
   
**fastqc.sh**: takes a single or list of fastq.gz and runs fastQC into an output folder called FASTQC_OUT
```
fastqc.sh "*gz"
```

**run_bbmerge**: Takes a single R1 fasta name and merges with R2.  Must be in casava style like \_R1_
```
run_bbmerge Sample_R1_001.fastq.gz
```

**cutadapt.sh**: trims off Illumina adapters and very low quality bases at the end (Q5) but can edit. Saves the output as Sample_cutadapt.txt so can be read in by multiqc.  Must be in \_R1_ format unless edited

```
cutadapt.sh Sample_R1_001.fastq.gz
```

## Assembly

**megahit.sh**: a VERY high memory default to take a list of R1 files with .R1. format.  Output is Basename.log, Basename.contigs.fa and Basename.
```
megahit_merge.sh *R1* #if you want all to be co-assembled
or
megahit_merge.sh Sample.R1.fastq.gz

```

**megahit_merge.sh**: very much like megahit but also expects a merge read file with the same name as .R1. but with .merge.

```
megahit_merge.sh *R1* if you want all to be co-assembled
or
megahit_merge.sh Sample.R1.fastq.gz
```

## Mapping reads

**bowtie2_index.sh**: runs bowtie2-build
```bowtie2_index.sh <input_fasta> <out_index_name>```

**bowtie2_map_sort.sh**: calls bowtie2 against a target genome, filters sam for mapped reads only, sorts and indexes.  Removes intermediate mapping files
```bowtie2_map_sort.sh <Read1> <Read2> <Output_Name> <Reference_index>```


## Combined pipeline for read processing, megahit, and mapping

**cut_assemble_map.sh**: does adapter trimming and quality trimming at Q20 with min len 31 and then assembles those reads with megahit and then maps them back to the assembly. Metrics are saved from the cutadapt and bowtie outputs.  Get trimmed fastq, contigs, and a sorted bam of ONLY MAPPED READS. Give it read 1.

```cut_assemble_map.sh <_R1_something.fastq.gz```

## Annotation
