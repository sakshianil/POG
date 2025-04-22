#!/bin/bash

# === CONFIGURATION ===
INPUT_DIR="data/proteomes"
OUTPUT_DIR="results/wideeffhunter"
WIDEEFFHUNTER_DIR="tools/WideEffHunter_1.0"
WIDEEFFHUNTER_SCRIPT="${WIDEEFFHUNTER_DIR}/WideEffHunter.sh"

# === ENSURE OUTPUT DIRECTORY EXISTS ===
mkdir -p "$OUTPUT_DIR"

# === LOOP THROUGH ALL PROTEOME FILES ===
for fasta in "$INPUT_DIR"/*.pep.all.fa; do
    species=$(basename "$fasta" .pep.all.fa)
    outdir="${OUTPUT_DIR}/${species}"

    echo " Running WideEffHunter for: $species"

    # Navigate to WideEffHunter script directory
    cd "$WIDEEFFHUNTER_DIR" || { echo " Cannot access $WIDEEFFHUNTER_DIR"; exit 1; }

    # Clean previous intermediate files
    rm -rf tmp Effectors input.fasta res_blast.txt *result.txt 2>/dev/null

    # Copy current proteome file to input.fasta
    cp "$PWD/../../$fasta" input.fasta

    # Run WideEffHunter
    bash "$WIDEEFFHUNTER_SCRIPT" ./input.fasta

    # Save results to designated output
    mkdir -p "$outdir"
    [ -d tmp ] && mv tmp "$outdir/"
    [ -d Effectors ] && mv Effectors "$outdir/"
    [ -f res_blast.txt ] && mv res_blast.txt "$outdir/"
    mv *result.txt "$outdir/" 2>/dev/null

    # Clean up input
    rm -f input.fasta

    echo " Completed WideEffHunter for: $species"
done

echo " All proteomes processed by WideEffHunter."

