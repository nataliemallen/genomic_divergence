#!/bin/bash
#SBATCH --job-name=mito_downstream
#SBATCH -A dewoody
#SBATCH -t 1-00:00:00
#SBATCH -p highmem
#SBATCH -n 96
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

# downstream mito trait analyses (STANDARD variant). Run this AFTER both PMM
# jobs (birds + mammals) have finished and written their scaled_df + tree
# objects + posteriors under /scratch/gautschi/allen715/GD_mito/standard/results/
### order matters: make_tree_objects rebuilds the tree1_objects that
# variance_partitioning needs. pPCA and the summary figures are independent

module load r/4.5.2
module load anaconda 
conda activate phylo

### echo "=== 1/4 rebuild tree1 objects ==="
# rscript mito_make_tree_objects.R
### echo "=== 2/4 variance partitioning ==="
# rscript mito_variance_paritioning.R

echo "=== 3/4 phylogenetic PCA (standard) ==="
Rscript mito_analyze_ppca.R

### echo "=== 4/4 summary tables + figures ==="
# rscript mito_analyze_all.R

echo "=== downstream complete ==="
