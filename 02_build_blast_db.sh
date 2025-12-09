ROOT="/Volumes/Public/SNP_array_data/horse_data"
REF="${ROOT}/ref/reference_genomes"

# Use the symlink that points to your local disk
DB_DIR="${ROOT}/blastdb"

FA="${REF}/EquCab3.0_Ensembl_toplevel.fa"

mkdir -p "${DB_DIR}"

makeblastdb \
  -in "${FA}" \
  -dbtype nucl \
  -parse_seqids \
  -out "${DB_DIR}/EquCab3_Ensembl"
