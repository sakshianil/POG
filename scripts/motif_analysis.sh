#!/bin/bash

# === Required paths ===
UPSTREAM_DIR="results/fasta/upstream_sequences"
PROTEIN_DIR="results/fasta/protein_sequences"
ALIGNMENT_DIR="results/alignments"
MOTIF_OUT="results/motif_analysis"
JASPAR_UP="resources/motif_databases/JASPAR2022_CORE_redundant_v2.meme"
ELM_PROT="resources/motif_databases/elm2024.meme"
GOMO_DB="resources/motif_databases/gomo_db"

mkdir -p "$MOTIF_OUT" "$ALIGNMENT_DIR"

# === MAFFT & Jalview alignments ===
for file in "$PROTEIN_DIR"/*.fasta; do
  cluster=$(basename "$file" .fasta)
  aln_file="$ALIGNMENT_DIR/${cluster}_aln.fasta"
  echo " Aligning $cluster with MAFFT"
  mafft --auto "$file" > "$aln_file"
	desktop_file="${aln_file%.fasta}.jv"
	echo "#MAFFT Alignment: $cluster" > "$desktop_file"
	echo "$aln_file" >> "$desktop_file"
  echo " Saved for Jalview: $desktop_file"
done

# === Motif Analysis ===
for up_file in "$UPSTREAM_DIR"/*.fasta; do
  cluster=$(basename "$up_file" .fasta)
  prot_file="$PROTEIN_DIR/${cluster}.fasta"
  cluster_dir="$MOTIF_OUT/$cluster"
  mkdir -p "$cluster_dir"

  # --- MEME Upstream ---
  meme "$up_file" -oc "$cluster_dir/MEME_oops_upstream" -mod oops -dna -nmotifs 10 -minw 6 -maxw 17 -evt 0.05 -markov_order 1
  meme "$up_file" -oc "$cluster_dir/MEME_zoops_upstream" -mod zoops -dna -nmotifs 10 -minw 6 -maxw 17 -markov_order 1

  # --- TOMTOM (DNA motifs) ---
  tomtom -oc "$cluster_dir/TOMTOM_oops_upstream" -thresh 0.05 -norc "$cluster_dir/MEME_oops_upstream/meme.xml" "$JASPAR_UP"
  tomtom -oc "$cluster_dir/TOMTOM_zoops_upstream" -thresh 0.05 -norc "$cluster_dir/MEME_zoops_upstream/meme.xml" "$JASPAR_UP"

  # --- AMA & GOMO ---
  for gomo_file in "$GOMO_DB"/*_1000_199.na; do
    gomo_base=$(basename "$gomo_file" .na)
    ama_out="$cluster_dir/ama_$gomo_base"
    gomo_out="$cluster_dir/gomo_$gomo_base"
    mkdir -p "$ama_out" "$gomo_out"

    ama "$cluster_dir/MEME_oops_upstream/meme.xml" "$gomo_file" "$GOMO_DB/${gomo_base}.na.bfile" --oc "$ama_out"

    gomo --verbosity 1 --oc "$gomo_out" --t 0.01 --gs --shuffle_scores 1000 --dag "$GOMO_DB/go.dag" \
      --motifs "$cluster_dir/MEME_oops_upstream/meme.xml" "$GOMO_DB/${gomo_base}.na.csv" "$ama_out/ama.xml"
  done

  # --- MEME Protein Motifs ---
  if [ -f "$prot_file" ]; then
    meme "$prot_file" -oc "$cluster_dir/MEME_oops_protein" -mod oops -nmotifs 10 -minw 6 -maxw 15 -evt 0.05 -markov_order 1
    meme "$prot_file" -oc "$cluster_dir/MEME_zoops_protein" -mod zoops -nmotifs 10 -minw 6 -maxw 15 -markov_order 1

    tomtom -oc "$cluster_dir/TOMTOM_oops_protein" -thresh 0.05 -norc "$cluster_dir/MEME_oops_protein/meme.xml" "$ELM_PROT"
    tomtom -oc "$cluster_dir/TOMTOM_zoops_protein" -thresh 0.05 -norc "$cluster_dir/MEME_zoops_protein/meme.xml" "$ELM_PROT"
  fi

echo " Completed motif analysis for $cluster"
done

