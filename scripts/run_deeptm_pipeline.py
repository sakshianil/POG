import os
import subprocess
from pathlib import Path
import pandas as pd
import argparse
import sys

# === PARSE ARGUMENTS ===
parser = argparse.ArgumentParser(description="Run DeepTM on secreted and non-secreted FASTA files.")
parser.add_argument('--signalp_dir', required=True, help='Path to SignalP output directory')
parser.add_argument('--deeptm_dir', required=True, help='Output directory for DeepTM results')
parser.add_argument('--npz_converter', required=True, help='Path to Fasta_to_deeptm_npz.py script')
parser.add_argument('--deeptm_runner', required=True, help='Path to DeepTM run_deep_tm.py script')
parser.add_argument('--model_dir', required=True, help='Path to saved DeepTM model directory')
args = parser.parse_args()

# === CONFIGURATION ===
signalp_dir = Path(args.signalp_dir)
deeptm_dir = Path(args.deeptm_dir)
npz_converter = Path(args.npz_converter)
inferencer = Path(args.deeptm_runner)
model_dir = Path(args.model_dir)
python_exec = sys.executable

deeptm_dir.mkdir(exist_ok=True)
summary_records = []

# Process each genome folder in signalp_output
for genome_dir in signalp_dir.glob("*/"):
    genome_name = genome_dir.name
    print(f"\n🔍 Processing genome: {genome_name}")

   # for category in ["secreted"]:  # Run only for secreted for now
    for category in ["secreted", "non_secreted"]: # Run for secreted and non_secreted   
        fasta_file = genome_dir / f"all_{category}.fasta"
        if not fasta_file.exists():
            print(f"⚠️  {category.capitalize()} file missing for {genome_name}, skipping.")
            continue

        output_subdir = deeptm_dir / genome_name / category
        output_subdir.mkdir(parents=True, exist_ok=True)

        npz_path = output_subdir / "input.npz"
        csv_path = output_subdir / "predictions.csv"
        npy_path = output_subdir / "raw_probs.npy"

        # Skip if already done
        if csv_path.exists() and npy_path.exists():
            print(f"✅ {category.capitalize()} already processed for {genome_name}, skipping.")
            continue

        print(f"➡️  Converting FASTA to .npz for {category}...")
        subprocess.run([
            python_exec, str(npz_converter),
            "--fasta", str(fasta_file),
            "--out", str(npz_path)
        ], check=True)

        print(f"➡️  Running DeepTM inference for {category}...")
        env = os.environ.copy()
        env["PYTHONPATH"] = str(inferencer.parent.resolve())
        subprocess.run([
            python_exec, str(inferencer.name),
            "--npz", str(npz_path.resolve()),
            "--model_dir", str(model_dir.resolve()),
            "--out_dir", str(output_subdir.resolve())
        ], check=True,cwd=str(inferencer.parent),env=env)

        # Add summary entry
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            counts = df['PredictedType'].value_counts().to_dict()
            row = {"Genome": genome_name, "Category": category}
            row.update(counts)
            summary_records.append(row)

# Save final summary
summary_df = pd.DataFrame(summary_records)
summary_df.fillna(0, inplace=True)
summary_path = deeptm_dir / "summary_deeptm.csv"
summary_df.to_csv(summary_path, index=False)
print(f"\n✅ DeepTM summary saved to {summary_path}")

