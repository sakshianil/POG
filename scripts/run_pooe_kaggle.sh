#!/bin/bash

# === CONFIGURATION ===
NOTEBOOK_PATH="scripts/pooe-lite.ipynb"
OUTPUT_DIR="results/pooe_lite_outputs"
OUTPUT_NOTEBOOK="${OUTPUT_DIR}/pooe-lite-executed.ipynb"

# === CREATE OUTPUT FOLDER ===
mkdir -p "$OUTPUT_DIR"

# === CHECK FOR PAPERMILL INSTALLATION ===
if ! command -v papermill &> /dev/null; then
    echo " Papermill not found. Please install it using 'pip install papermill'."
    exit 1
fi

# === EXECUTE NOTEBOOK ===
echo "Running POOE Lite notebook with papermill..."
papermill "$NOTEBOOK_PATH" "$OUTPUT_NOTEBOOK"

if [ $? -eq 0 ]; then
    echo "Successfully ran pooe-lite notebook. Output saved to: $OUTPUT_NOTEBOOK"
else
    echo " Notebook execution failed."
    exit 1
fi

