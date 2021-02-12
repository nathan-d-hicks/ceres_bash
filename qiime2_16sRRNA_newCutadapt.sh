#!/bin/bash

#SBATCH --time=47:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=6   # 1 processor core(s) per node X 2 threads per core
#SBATCH --mem=30G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="qiime_assign"
#SBATCH --mail-user=nathan.hicks@usda.gov   # email address
#SBATCH --mail-type=END

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#Define the input path for the data
INPATH='/project/epicon/nhicks/syncom2/raw'
OUTPATH='/project/epicon/nhicks/syncom2/qiime2/modify_cutadapt'
# need to set to somewhere you can write a lot of data!
TMPDIR='/project/epicon/nhicks/syncom2/qiime2/'

#Define path to classifier
CLASSIFIER='/home/nathan.hicks/silva_database/silva-138-99-341F-785R-classifier.qza'

module load qiime2/2020.6

#Import data -- IMPORTANT NO UNDERSCORES IN SAMPLE
#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $INPATH --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path ${OUTPATH}/demux-paired.qza 

#summarize the stats of reads and quality
#qiime demux summarize --i-data $OUTPATH/demux-paired.qza --o-visualization $OUTPATH/demux.qzv

#Remove adapter pairs
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $OUTPATH/demux-paired.qza \
  --p-cores 6 \
  --p-front-f CCTACGGGNBGCASCAG \
  --p-front-r GACTACNVGGGTATCTAATCC \
  --p-overlap 3 \
  --p-times 2 \
  --o-trimmed-sequences $OUTPATH/demux-paired-trimmed.qza \
  --output-dir $OUTPATH/unspecified-trim \
  --verbose \
  --p-discard-untrimmed

#Get a summary after trimming
qiime demux summarize \
  --i-data $OUTPATH/demux-paired-trimmed.qza \
  --o-visualization $OUTPATH/demux-trimmed.qzv

## check you read length and quality using the .qzv files to determine where you want to trim your reads
## reads needs to be long enough to overlap when joining paired ends
## note: set --p-n-reads-learn based on your read length
## a million 100nt reads (or 100M total bases) is adequate to learn the error rates
## e.g. 100M/280 ~ 360000 --p-n-reads-learn 360000
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $OUTPATH/demux-paired-trimmed.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 200 \
  --p-n-threads 0 \
  --p-n-reads-learn 360000 \
  --o-table $OUTPATH/table.qza \
  --o-representative-sequences $OUTPATH/rep-seqs.qza \
  --o-denoising-stats $OUTPATH/denoising-stats.qza \
  --output-dir $OUTPATH/unspecified-dada2 \
  --verbose

qiime feature-table summarize \
  --i-table $OUTPATH/table.qza \
  --o-visualization $OUTPATH/table.qzv \
#  --m-sample-metadata-file $OUTPATH/sample-metadata-final_feature.tsv

qiime feature-table tabulate-seqs \
  --i-data $OUTPATH/rep-seqs.qza \
  --o-visualization $OUTPATH/rep-seqs.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $OUTPATH/rep-seqs.qza \
  --o-alignment $OUTPATH/aligned-rep-seqs.qza \
  --o-masked-alignment $OUTPATH/masked-aligned-rep-seqs.qza \
  --o-tree $OUTPATH/unrooted-tree.qza \
  --o-rooted-tree $OUTPATH/rooted-tree.qza

qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads $OUTPATH/rep-seqs.qza \
  --o-classification $OUTPATH/taxonomy.qza

## generate a summary
qiime metadata tabulate \
  --m-input-file $OUTPATH/taxonomy.qza \
  --o-visualization $OUTPATH/taxonomy.qzv

## exporting data
## exporting a feature table
qiime tools export \
  --input-path $OUTPATH/table.qza \
  --output-path $OUTPATH/exported-feature-table

## exporting a phylogenetic tree
qiime tools export \
  --input-path $OUTPATH/unrooted-tree.qza \
  --output-path $OUTPATH/exported-tree

## exporting a taxonomy table
qiime tools export \
  --input-path $OUTPATH/taxonomy.qza \
  --output-path $OUTPATH/exported-taxonomy
