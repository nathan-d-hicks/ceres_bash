#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12-12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=24  # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=100G
#SBATCH --job-name="kaiju"
#SBATCH --partition=priority-mem
#SBATCH --qos=pgec

module load kaiju

kaiju -t /reference/data/kaiju/2020-04-02/nodes.dmp -f /reference/data/kaiju/2020-04-02/nr/kaiju_db_nr.fmi -o $2 -z 24 \
    -i $1