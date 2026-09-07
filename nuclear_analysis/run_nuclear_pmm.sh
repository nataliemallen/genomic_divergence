#!/bin/bash
#SBATCH --job-name=nuclear_pmm
#SBATCH -A dewoody
#SBATCH -t 1-00:00:00
#SBATCH -p highmem
#SBATCH -n 128
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

# nuclear dyadic phylogenetic mixed model
# to run: sbatch run_nuclear_pmm.sh birds
# sbatch run_nuclear_pmm.sh mammals

module load r/4.5.2

export RUN_CLASS="${1:-birds}"
echo "=== nuclear pmm | class = ${RUN_CLASS} ==="
Rscript nuclear_pmm.R

