# Run this locally before sbatch
curl "https://docs.google.com/spreadsheets/d/1FUBfBHiFseFQ5iJKXaoXmbx3vbuBNreRsrCyLk2p9qg/export?format=csv&gid=0" \
  -L -o candidate_genes_extra.csv
echo "" >>candidate_genes_extra.csv
scp candidate_genes_extra.csv renaut@fir.alliancecan.ca:/project/def-rallard/COMMUN/raredisease_rnaseq/data/input/.
rm candidate_genes_extra.csv
