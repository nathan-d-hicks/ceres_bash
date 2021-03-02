# 16s Qiime Analysis with DADA2

**qiime2_16sRRNA.sh**: This is a pipeline to run Qiime from paired, demultiplexed reads through read merging, ASV calling, taxonimic identification and output of a tree, taxonomic assignment, and readcount table. Requires editing to change in input path, output path, temp directory, and classifier location.


**qiime2_16sRRNA_newCutadapt.sh**: This is a pipeline identical to the one above, except that cutadapt parameters have been modified to trim more aggressively including allowing multiple copies of the adapter in a row. This seems to deal with a common problem where you get multiple ASVs that are identical except one has an adapter stuck on the end of it. It seems like you sometimes might really get 2 adapters in a row which maybe is a PCR artifact?? It works.

**qiime_taxonomy.sh**: Runs only the naive bayesian classifier portion of the qiime2 pipeline.  Hard-coded to use the silva 138 99percent database in /project/epicon/bin/silva_database/  Output is an exported taxonomy table.
```
qiime_taxonomy.sh <reads.fasta> <output_name>
```



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

**split_interleaved.sh**: Splits an interleaved fastq file into a forward and reverse fastq.  Hardcoded to expect the file to end in \_filtered.fastq but this can be easily edited.

```
split_interleaved.sh Sample_interleaved_filtered.fastq
```

**run_sortmerna.sh**: Aligns paired reads against all of the sortmerna databases stored in /project/epicon/bin/.  This again is hard-coded to expect the input file to end in \_1.trimmed.fastq.gz but can be easily edited. Requires a LOT of memory for reasons I don't understand.  Outputs a file rRNA.fastq which is pairs where either read mapped to rRNA and then filtered.fastq which is pairs where neither aligned to rRNA.
```
run_sortmerna.sh Sample_1.trimmed.fastq.gz
```


## Assembly

**megahit.sh**: a VERY high memory default to take a list of R1 files with .R1. format.  Output is Basename.log, Basename.contigs.fa and Basename.
```
megahit_merge.sh *R1* #if you want all to be co-assembled
or
megahit_merge.sh Sample.R1.fastq.gz

```

**megahit_merge.sh**: very much like megahit but also expects a merge read file with the same name as .R1. and .R2. but with .merge.

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

**run_kallisto.sh**: An essentially useless script because it uses my conda environments and is hard coded to the genome I was interested in. But is a good template for how to call kallisto I guess.


## Combined pipeline for read processing, megahit, and mapping

**cut_assemble_map.sh**: does adapter trimming and quality trimming at Q20 with min len 31 and then assembles those reads with megahit and then maps them back to the assembly. Metrics are saved from the cutadapt and bowtie outputs.  Get trimmed fastq, contigs, and a sorted bam of ONLY MAPPED READS. Give it read 1.

```cut_assemble_map.sh <_R1_something.fastq.gz```

**raw_to_assembly_pipeline.sh**: Same as above but before assembly, it aligns the the trimmed reads to a hardmasked copy of the sorghum genome.  Then read pairs where neither mapped to sorghum are extracted and used for all subsequent steps.
```
raw_to_assembly_pipeline.sh <_R1_something.fastq.gz
```

## Annotation and Taxonomic Assignments


**run_prokka.sh**: Calls prokka on a metagenome with default parameters.
```
run_prokka.sh Metagenome.fasta
```

**run_prokka_sg.sh**: Calls prokka on not-a-metagenome. Good for single genomes or bins.  The output name is a little more sophisticated so that it appends \_prokka to the end of all of the text before the last period.
```
run_prokka_sg.sh Genome.fasta
```

**kaiju.sh**: Calls a pretty vanilla version of kaiju to do taxonomic classification based on a verions of the nr database that was downloaded on April 2 2020 by....somebody. You can find the location of the database by looking in the script.
```
kaiju.sh <contigs.fasta> <outname>
```
