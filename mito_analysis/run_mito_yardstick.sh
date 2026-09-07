#!/bin/bash
#SBATCH --job-name=mito_yardstick
#SBATCH -A dewoody
#SBATCH -t 08:00:00
#SBATCH -p cpu
#SBATCH -n 32
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

# mitochondrial yardstick: class-level + order-level. No trait analysis
# prereq: place Complete_divergence_with_mito_05-14-26.csv in
# scratch/gautschi/allen715/GD_models/final_dataset
# (the scripts read it from there, alongside the traits .xlsx)
# outputs are written under /scratch/gautschi/allen715/GD_models/ in
# mito_yardstick/ and mito_order_yardstick/ (separate from nuclear results)

module load r/4.5.2

echo "=== class-level mito yardstick ==="
Rscript mito_yardstick.R

echo "=== order-level mito yardstick ==="
# rscript mito_order_yardstick.R

echo "=== done ==="
