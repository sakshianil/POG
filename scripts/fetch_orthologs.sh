#!/bin/bash

# === CONFIGURATION ===
INPUT_FILE="data/final_counts_Pl_halstedii.txt"         # Input gene counts from Pl. halstedii transcriptome from project (PRJEB49134)
OUTPUT_FILE="results/orthologs/ortholog_summary.csv"    # Output would contain all orthologs found within SAR supergroup. Select >70% seqeunce identity only for further analsyis
TEMP_FILE="results/orthologs/temp_results.csv"

mkdir -p "$(dirname "$OUTPUT_FILE")"
echo "Gene_ID,Matched_Gene,Method_Link,Protein_ID,Species,Taxonomy_Level,Ortholog_Type,Source_perc_id,Source_perc_pos,Target_perc_id,Target_perc_pos" > "$OUTPUT_FILE"

# === START PROCESSING ===
tail -n +2 "$INPUT_FILE" | cut -f1 | while read gene; do
    echo "Processing: $gene"

    # Get ortholog data
    response=$(wget -q --header='Content-type:text/xml' "https://rest.ensembl.org/homology/id/${gene}?type=orthologues;sequence=protein" -O -)

    # If empty response
    if [[ -z "$response" ]]; then
        echo "$gene,NO_MATCH,NO_MATCH,NO_MATCH,NO_MATCH,NO_MATCH,NO_MATCH,NO_MATCH,NO_MATCH,NO_MATCH,NO_MATCH" >> "$OUTPUT_FILE"
        continue
    fi

    # Extract values from XML
    method_link=$(echo "$response" | xmllint --xpath 'string(//homologies/@method_link_type)' - 2>/dev/null)
    taxonomy_level=$(echo "$response" | xmllint --xpath 'string(//homologies/@taxonomy_level)' - 2>/dev/null)
    ortholog_type=$(echo "$response" | xmllint --xpath 'string(//homologies/@type)' - 2>/dev/null)

    query_gene=$(echo "$response" | xmllint --xpath 'string(//homologies/source/@id)' - 2>/dev/null)
    source_perc_id=$(echo "$response" | xmllint --xpath 'string(//homologies/source/@perc_id)' - 2>/dev/null)
    source_perc_pos=$(echo "$response" | xmllint --xpath 'string(//homologies/source/@perc_pos)' - 2>/dev/null)

    target_gene=$(echo "$response" | xmllint --xpath 'string(//homologies/target/@id)' - 2>/dev/null)
    protein_id=$(echo "$response" | xmllint --xpath 'string(//homologies/target/@protein_id)' - 2>/dev/null)
    species=$(echo "$response" | xmllint --xpath 'string(//homologies/target/@species)' - 2>/dev/null)
    target_perc_id=$(echo "$response" | xmllint --xpath 'string(//homologies/target/@perc_id)' - 2>/dev/null)
    target_perc_pos=$(echo "$response" | xmllint --xpath 'string(//homologies/target/@perc_pos)' - 2>/dev/null)

    # Fallback values
    query_gene="${query_gene:-$gene}"
    target_gene="${target_gene:-NO_MATCH}"
    method_link="${method_link:-NO_MATCH}"
    protein_id="${protein_id:-NO_MATCH}"
    species="${species:-NO_MATCH}"
    taxonomy_level="${taxonomy_level:-NO_MATCH}"
    ortholog_type="${ortholog_type:-NO_MATCH}"
    source_perc_id="${source_perc_id:-NO_MATCH}"
    source_perc_pos="${source_perc_pos:-NO_MATCH}"
    target_perc_id="${target_perc_id:-NO_MATCH}"
    target_perc_pos="${target_perc_pos:-NO_MATCH}"

    # Write to file
    echo "$query_gene,$target_gene,$method_link,$protein_id,$species,$taxonomy_level,$ortholog_type,$source_perc_id,$source_perc_pos,$target_perc_id,$target_perc_pos" >> "$OUTPUT_FILE"

    sleep 1  # Avoid rate limiting
done

echo " Ortholog fetch completed. Output saved to: $OUTPUT_FILE"

