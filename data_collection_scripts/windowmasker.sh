!/bin/bash
#SBATCH --job-name=windowmasker
#SBATCH -A dewoody
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -q normal
#SBATCH -p cpu
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load samtools

export PATH=$PATH:/home/allen715/bin/winmasker/ncbi_cxx--25_2_0/GCC850-ReleaseMT64/bin/

CSV="batch2_genomes.csv"
GENOME_DIR="batch2_final_genomes"
OUTDIR="windowmasker_files"
RESULTS="repeat_content_summary.csv"

mkdir -p "$OUTDIR"

# results CSV
echo "accession,repeat_content_percent" > "$RESULTS"

# read accessions from CSV csv
tail -n +2 "$CSV" | while IFS=, read -r accession _; do
  accession=$(echo "$accession" | xargs)
  genome="${GENOME_DIR}/${accession}_final.fna"

wm_counts="${OUTDIR}/${accession}.wm.counts"
  masked="${OUTDIR}/${accession}_masked.fna"

  # skip if output already exists
  if [[ -f "$masked" || -f "$wm_counts" ]]; then
    echo "skipping $accession - output exists"
    continue
  fi

  if [[ ! -f "$genome" ]]; then
    echo "missing genome file $genome" >&2
    continue
  fi

  # generate counts
  windowmasker -mk_counts -in "$genome" -out "$wm_counts"

  # mask genome
  windowmasker -ustat "$wm_counts" -in "$genome" -out "$masked" -outfmt fasta -dust true

  # calculate repeat content
  repeat_pct=$(grep -v ">" "$masked" | awk '
  {
    seq = toupper($0);
    total += length(seq);
    masked += gsub(/[a-z]/, "", $0);
  }
  END { if (total > 0) printf "%.2f", 100 * masked / total; else print "0.00" }')

  echo "$accession,$repeat_pct" >> "$RESULTS"
  echo "$accession done — repeat content: $repeat_pct%"
done
