import os
import subprocess
import logging
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
import pandas as pd
import argparse

# Add at the top after imports
parser = argparse.ArgumentParser(description="Run SignalP6 on a directory of proteomes.")
parser.add_argument('--input_dir', required=True, help='Directory with .pep.all.fa proteome FASTA files')
parser.add_argument('--output_dir', required=True, help='Directory to write SignalP outputs')
parser.add_argument('--min_seq_length', type=int, default=50, help='Minimum sequence length to include')
args = parser.parse_args()

# Replace hardcoded paths with:
input_dir = args.input_dir
output_root = Path(args.output_dir)
min_seq_length = args.min_seq_length


# SETUP LOGGING
output_root.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=output_root / 'signalp_pipeline.log',
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logging.info("Started SignalP pipeline.")

# CHECK SIGNALP INSTALLATION
def check_signalp_installed():
    try:
        subprocess.run(["signalp6", "--version"], check=True, stdout=subprocess.DEVNULL)
        logging.info("SignalP 6.0 detected.")
    except Exception as e:
        logging.error("SignalP 6.0 is not installed or not in PATH.")
        raise RuntimeError("SignalP 6.0 is not installed or not in PATH.") from e

# PARSE SIGNALP OUTPUT TO GET SECRETED IDS
def get_secreted_ids(prediction_file):
    secreted_ids = set()
    with open(prediction_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2 and parts[1].strip() != "OTHER":
                secreted_ids.add(parts[0])
    return secreted_ids

# PROCESS EACH FILE
def process_proteome(fasta_file):
    basename = fasta_file.stem.replace(".pep.all", "")
    output_dir = output_root / basename
    output_dir.mkdir(parents=True, exist_ok=True)

    prediction_file = output_dir / "prediction_results.txt"
    secreted_out = output_dir / "all_secreted.fasta"
    non_secreted_out = output_dir / "all_non_secreted.fasta"

    # Skip if already processed
    if prediction_file.exists() and secreted_out.exists() and non_secreted_out.exists():
        logging.info(f"Skipping {basename} — already processed.")
        return (basename, -1, -1, "SKIPPED")

    # Filter input sequences
    filtered_fasta = output_dir / "filtered_input.fasta"
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    long_sequences = [seq for seq in sequences if len(seq.seq) >= min_seq_length]
    SeqIO.write(long_sequences, filtered_fasta, "fasta")
    logging.info(f"{basename}: Filtered {len(sequences) - len(long_sequences)} short sequences.")

    # Run SignalP
    try:
        subprocess.run([
            "signalp6",
            "--fastafile", str(filtered_fasta),
            "--organism", "other",
            "--output_dir", str(output_dir),
            "--format", "txt",
            "--mode", "fast"
        ], check=True)
        logging.info(f"{basename}: SignalP prediction complete.")
    except subprocess.CalledProcessError as e:
        logging.error(f"{basename}: SignalP failed - {e}")
        return (basename, -1, -1, "ERROR")

    # Parse prediction results
    secreted_ids = get_secreted_ids(prediction_file)
    all_ids = set(rec.id for rec in long_sequences)

    # Write FASTA files
    secreted_records = []
    non_secreted_records = []
    for rec in long_sequences:
        if rec.id in secreted_ids:
            secreted_records.append(rec)
        else:
            non_secreted_records.append(rec)

    SeqIO.write(secreted_records, secreted_out, "fasta")
    SeqIO.write(non_secreted_records, non_secreted_out, "fasta")

    logging.info(f"{basename}: Wrote {len(secreted_records)} secreted, {len(non_secreted_records)} non-secreted sequences.")
    return (basename, len(secreted_records), len(non_secreted_records), "DONE")

# MAIN
def main():
    check_signalp_installed()
    results = []

    for fasta in Path(input_dir).glob("*.pep.all.fa"):
        result = process_proteome(fasta)
        results.append(result)

    # Generate summary
    summary_df = pd.DataFrame(results, columns=["Proteome", "Secreted", "Non_Secreted", "Status"])
    summary_path = output_root / "summary_signalp.csv"
    summary_df.to_csv(summary_path, index=False)
    logging.info(f"Summary written to {summary_path}")

if __name__ == "__main__":
    main()
