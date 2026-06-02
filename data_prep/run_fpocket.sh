#!/bin/bash

# Script to run fpocket prediction on protein structures
# install fpocket with conda first: conda create -n fpocket_env -c bioconda fpocket

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate fpocket_env

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Parent directory (the repo root)
REPO_DIR="$(dirname "$SCRIPT_DIR")"
# Data directory
DATA_DIR="$REPO_DIR/data"
# Input directory with PDB frames from bioemulations
RESULTS_BASE_DIR="$DATA_DIR/bioemu_results"
# Output directory for fpocket predictions
OUTPUT_BASE_DIR="$DATA_DIR/fpocket_preds"
# Number of parallel processes (adjust based on your system)
NUM_PROCESSES=32

# Check if fpocket executable exists
if ! command -v fpocket &> /dev/null; then
    echo "Error: fpocket executable not found in the current environment."
    echo "Please make sure fpocket is installed and the conda environment is activated."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_BASE_DIR"

# Counter for tracking progress
total=0
processed=0

# Export variables for parallel execution
export RESULTS_BASE_DIR
export OUTPUT_BASE_DIR

# Function to process a single protein directory
process_protein() {
    protein_dir=$1
    protein_id=$(basename "$protein_dir")
    frames=("$protein_dir"/frame_*.pdb)

    # If no matching frame files exist, skip this directory
    if [ ! -f "${frames[0]}" ]; then
        return
    fi

    echo "Processing protein: $protein_id"
    protein_output_dir="$OUTPUT_BASE_DIR/$protein_id"
    mkdir -p "$protein_output_dir"

    # Run fpocket on each frame
    for frame in "${frames[@]}"; do
        fpocket -f "$frame"
        # fpocket creates output in the same directory as the input file.
        # Let's move the results to our output directory.
        frame_base=$(basename "$frame" .pdb)
        output_dir_name="${frame_base}_out"
        mv "$protein_dir/$output_dir_name" "$protein_output_dir/"
    done
    echo "  Completed: $protein_id"
}

export -f process_protein

# Find all protein directories
protein_dirs=()
for d in "$RESULTS_BASE_DIR"/*; do
    if [ -d "$d" ]; then
        protein_dirs+=("$d")
    fi
done
total=${#protein_dirs[@]}

echo "Found $total proteins to process."

# Run in parallel
printf "%s\n" "${protein_dirs[@]}" | xargs -I {} -P "$NUM_PROCESSES" bash -c 'process_protein "{}"'

processed=$(find "$OUTPUT_BASE_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)

echo ""
echo "========================================="
echo "fpocket prediction complete!"
echo "Processed: $processed/$total proteins"
echo "Results saved to: $OUTPUT_BASE_DIR"
echo "========================================="
