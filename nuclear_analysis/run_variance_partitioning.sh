#!/bin/bash
#SBATCH --job-name=bird_phylo
#SBATCH -A dewoody
#SBATCH -t 1-00:00:00 
#SBATCH -p cpu
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load r/4.5.2

Rscript variance_partitioning_consensus.R
