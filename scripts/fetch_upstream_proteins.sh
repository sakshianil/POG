#!/bin/bash

# === INPUT FILE ===
INPUT_FILE="data/selected_orthologs.txt"  # Tab-delimited: Cluster_ID<TAB>Gene_ID selected from orthologs_summary.csv

# === OUTPUT DIRECTORIES ===
BASE_DIR="results/fasta"
PROTEIN_DIR="$BASE_DIR/protein_sequences"
UPSTREAM_DIR="$BASE_DIR/upstream_sequences"

mkdir -p "$PROTEIN_DIR" "$UPSTREAM_DIR"

# === FUNCTION TO DETERMINE SPECIES NAME FROM GENE ID ===
get_species_name() {
    local gene_id="$1"
    case "$gene_id" in
        Hpa*) echo "hyaloperonospora_arabidopsidis" ;;
        PYU1*) echo "pythium_ultimum" ;;
        PHYSODRAFT*) echo "phytophthora_sojae" ;;
        PITG*) echo "phytophthora_infestans" ;;
        *) echo "plasmopara_halstedii_gca_900000015" ;;
    esac
}

# === MAIN LOOP ===
awk -F'\t' 'NR > 1 {print $1, $2}' "$INPUT_FILE" | sort -u | while read cluster_id gene_id; do
    echo " Cluster: $cluster_id | Gene: $gene_id"

    # Output FASTA headers
    fasta_header=">${gene_id}"

    protein_fasta="${PROTEIN_DIR}/Cluster_${cluster_id}_proteins.fasta"
    upstream_fasta="${UPSTREAM_DIR}/Cluster_${cluster_id}_upstream.fasta"

    ### FETCH PROTEIN SEQUENCE ###
    protein_seq=$(wget -q -O - "https://rest.ensembl.org/sequence/id/${gene_id}?type=protein&content-type=text/plain")

    if [[ -n "$protein_seq" ]]; then
        echo "$fasta_header" >> "$protein_fasta"
        echo "$protein_seq" >> "$protein_fasta"
    else
        echo " No protein found for $gene_id"
    fi

    ### FETCH GENE INFO ###
    gene_info=$(wget -q -O - "https://rest.ensembl.org/overlap/id/${gene_id}?feature=gene&content-type=application/json")
    strand_gene=$(echo "$gene_info" | jq -r '.[0].strand')
    gene_start=$(echo "$gene_info" | jq -r '.[0].start')
    gene_end=$(echo "$gene_info" | jq -r '.[0].end')
    seq_region_name=$(echo "$gene_info" | jq -r '.[0].seq_region_name')

    # Calculate upstream coordinates
    if [[ "$strand_gene" == "1" ]]; then
        upstream_start=$((gene_start - 1000))
        [[ $upstream_start -lt 1 ]] && upstream_start=1
        upstream_end=$((gene_start - 1))
    else
        upstream_start=$((gene_end + 1))
        upstream_end=$((gene_end + 1000))
    fi

    # Get species
    species_name=$(get_species_name "$gene_id")

    ### FETCH UPSTREAM SEQUENCE ###
    upstream_seq=$(wget -q -O - \
        "https://rest.ensembl.org/sequence/region/${species_name}/${seq_region_name}:${upstream_start}:${upstream_end}:${strand_gene}?content-type=text/plain")

    if [[ -n "$upstream_seq" ]]; then
        echo "$fasta_header" >> "$upstream_fasta"
        echo "$upstream_seq" >> "$upstream_fasta"
    else
        echo " No upstream sequence found for $gene_id"
    fi

    sleep 1  # Be kind to the API
done

echo " Upstream and protein sequence fetch complete."

