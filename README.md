# Nextflow Pipeline for Oomycete Effector Analysis

## Overview
This pipeline performs comprehensive analysis on five oomycete proteomes to identify secreted effectors and investigate their regulatory architecture. It integrates tools for secretion prediction, effector classification, ortholog identification, upstream sequence extraction, motif discovery, and expression integration. The nextflow pipeline (**main.nf**) contains the further information on step-by-step scripts usage within pipeline. 

## Workflow Summary

### 1. **Download Proteomes**
- Source: Ensembl Protists Release 59
- Using ftp for all genomes (For example:  "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-60/protists/fasta/protists_stramenopiles1_collection/plasmopara_halstedii_gca_900000015/pep/Plasmopara_halstedii_gca_900000015.Plasmopara_halstedii_genome.pep.all.fa.gz"). Download five proteomes in data/proteomes
- Species:
  - *Plasmopara halstedii* (GCA_900000015.1)
  - *Phytophthora sojae* (GCA_000149755.1)
  - *Phytophthora infestans* (GCA_000142945.1)
  - *Hyaloperonospora arabidopsidis* (GCA_000173235.2)
  - *Pythium ultimum* (GCA_000143045.1)

### 2. **Secretion Signal Prediction**
- Tool: SignalP 6.0 
- Script used:scripts/signalp_pipeline.py 
- Output: secreted vs. non-secreted FASTA in directory results/

### 3. **Transmembrane Domain Prediction**
- Tool: DeepTM
- Script used:scripts/deeptm_pipeline.py
- Classification: SP+Glob, SP+TM, TM, Glob

### 4. **Effector Prediction**
- Tools:
  - POOE (via Kaggle notebook in repo: https://github.com/sakshianil/POOE_lite)
  - WideEffHunter v1.0 ( Script used: scripts/run_wideeffhunter.sh)
- Output: Effectors grouped as motif/domain/PHI hits/cysteine-rich

### 5. **Ortholog Identification**
- Source: Ensembl REST API
- Type: One-to-one orthologs only 
- Script used: scripts/fetch_orthologs.sh
- Validation: MAFFT + Jalview + InterProScan domains

### 6. **Upstream and Protein Sequence Retrieval**
- Tool: Ensembl REST API
- Script used: scripts/fetch_upstream_protein.sh 
- Region: 1000 bp upstream from TSS

### 7. **Motif Discovery and Annotation**
- Tools: MEME Suite (MEME, FIMO, MAST, GOMO, AMA, TOMTOM)
- Script used: scripts/motif_analysis.sh
- Databases stored within resources folder (define path to all databases within script as well):
  - JASPAR 2024
  - ELM 2024
  - Custom GOMO dbs

### 8. **Expression Data Integration**
- Dataset: *Pl. halstedii* transcriptome (PRJEB49134)
- Tool: DESeq2 (pipeline adapted from https://github.com/sakshianil/Transcriptional_regulation_oomycetes)
- Use: Match conserved orthologs to time-series profiles

### 9. **Visualization**
- Tools:
  - R (ggplot2, pheatmap)
  - motifStack (for motif phylogeny and logos)
- Output: Clustered heatmaps and motif enrichment maps

## Acknowledgments
- POOE: https://github.com/zzdlabzm/POOE
- Structural layout influenced by OpenAI support for reproducible pipelines.

## GitHub Repository
https://github.com/sakshianil/POOE_lite
https://github.com/sakshianil/Transcriptional_regulation_oomycetes

## License
MIT License

## How to Run
– The Nextflow `main.nf` script and `config` file 

