#!/bin/bash
#SBATCH --job-name=mito_align
#SBATCH -A dewoody
#SBATCH -t 6-00:00:00 
#SBATCH -p cpu
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module --force purge
module load biocontainers
module load mafft/7.490
module load anaconda

export NCBI_EMAIL="" 
  
python download_mitogenomes.py \
    --species species_list.csv \
    --pairs   pairs.csv \
    --outdir  results
    
python align_and_k2p.py --coverage results/pair_coverage.csv \
    --fasta results/fasta --cds results/cds --outdir results --threads 4 \
    --mafft-cmd "singularity exec --bind $PWD /apps/biocontainers/images/quay.io_biocontainers_mafft:7.490--h779adbc_0.sif mafft"

python merge_mito.py \
    --divergence Complete_divergence_05-14-26.csv \
    --compare    results/mito_k2p_compare.csv \
    --out        Complete_divergence_with_mito_05-14-26.csv