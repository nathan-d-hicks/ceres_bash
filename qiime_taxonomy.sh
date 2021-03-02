#!/bin/bash

#SBATCH --time=47:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # 1 processor core(s) per node X 2 threads per core
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="qiime_assign"
#SBATCH --mail-user=nathan.hicks@usda.gov   # email address
#SBATCH --mail-type=END

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load qiime2/2020.6

inputFasta=$1
outputName=$2
outputSeqs=${outputName}.qza
outputTax=${outputName}.taxonomy.qza
outputExport=${outputName}.taxonomy-export


qiime tools import --type FeatureData[Sequence] --input-path $1 --input-format DNAFASTAFormat --output-path $outputName

qiime feature-classifier classify-sklearn --i-classifier /project/epicon/bin/silva_database/silva-138-99-341F-785R-classifier.qza --i-reads $outputSeqs --o-classification $outputTax --p-reads-per-batch 2000

qiime tools export --output-path $outputExport --input-path $outputTax
