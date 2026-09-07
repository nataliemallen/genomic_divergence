#!/bin/bash
#SBATCH --job-name=mito_pmm
#SBATCH -A dewoody
#SBATCH -t 1-00:00:00
#SBATCH -p highmem
#SBATCH -n 128
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

# how to run
# sbatch run_mito_pmm.sh birds
# sbatch run_mito_pmm.sh mammals
# outputs -> /scratch/gautschi/allen715/GD_mito/standard/results
# resubmit will resume

module load r/4.5.2

export RUN_CLASS="${1:-birds}"
echo "=== mito PMM (standard) | class = ${RUN_CLASS} ==="
Rscript mito_pmm.R

