#!/bin/bash

# Define input and output paths
INPUT_DIR="data/proteomes"
SCRIPT="scripts/run_signalp_pipeline.py"
OUTPUT_DIR="results/signalp_output"

# Create output directory if not exists
mkdir -p "$OUTPUT_DIR"

# Log start
echo " Running SignalP6 on all proteomes from $INPUT_DIR..."
echo "Using script: $SCRIPT"

# Run the script
python3 "$SCRIPT" \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --min_seq_length 50  # You can adjust this if needed

# Check exit code
if [ $? -eq 0 ]; then
    echo "SignalP6 completed successfully. Output in $OUTPUT_DIR"
else
    echo "SignalP6 encountered an error."
fi

