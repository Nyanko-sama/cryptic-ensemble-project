#!/bin/bash

# Script to run p2rank prediction on protein structures
# Requires Java to be installed and p2rank to be downloaded and extracted (run get_p2rank.sh first)
# On metacentrum: module add openjdk/17

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Parent directory (the repo root)
REPO_DIR="$(dirname "$SCRIPT_DIR")"
# Data directory
DATA_DIR="$REPO_DIR/data"
# P2rank executable
P2RANK_BIN="$REPO_DIR/p2rank/prank"
# Results directory
RESULTS_BASE_DIR="$DATA_DIR/bioemu_results"
# Output directory
OUTPUT_BASE_DIR="$DATA_DIR/p2rank_preds"
# Number of parallel processes (adjust based on your system)
NUM_PROCESSES=4

# Check if p2rank executable exists
if [ ! -f "$P2RANK_BIN" ]; then
    echo "Error: p2rank executable not found at $P2RANK_BIN"
    echo "Please run get_p2rank.sh first to download and extract p2rank"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_BASE_DIR"

# Counter for tracking progress
total=0
processed=0

# Iterate over each protein directory
for protein_dir in "$RESULTS_BASE_DIR"/*; do
    if [ ! -d "$protein_dir" ]; then
        continue
    fi

    protein_id=$(basename "$protein_dir")
    frames=("$protein_dir"/frame_*.pdb)

    # If no matching frame files exist, skip this directory
    if [ ! -f "${frames[0]}" ]; then
        continue
    fi

    total=$((total + 1))
    protein_output_dir="$OUTPUT_BASE_DIR/$protein_id"
    mkdir -p "$protein_output_dir"

    dataset_file="$protein_output_dir/dataset.ds"

    echo "Processing protein: $protein_id"
    echo "  Writing dataset file with ${#frames[@]} frames..."

    # Create a dataset file listing one PDB path per line
    : > "$dataset_file"
    for frame in "${frames[@]}"; do
        echo "$frame" >> "$dataset_file"
    done

    echo "  Running p2rank prediction with alphafold mode..."
    "$P2RANK_BIN" predict -o "$protein_output_dir" -c alphafold -t "$NUM_PROCESSES" "$dataset_file"
    if [ $? -eq 0 ]; then
        ((processed++))
        echo "  Completed"
    else
        echo "  Failed"
    fi

done

echo ""
echo "========================================="
echo "P2rank prediction complete!"
echo "Processed: $processed/$total proteins"
echo "Results saved to: $OUTPUT_BASE_DIR"
echo "========================================="
