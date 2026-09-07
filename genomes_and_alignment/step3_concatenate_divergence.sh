echo "pair,raw,k2p,k3p,jc" > batch1_divergence_update.csv

for file in *_distance.csv; do
    tail -n +2 "$file" >>batch1_divergence_update.csv
done
