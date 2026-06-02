#!/bin/bash

# Script to run p2rank fpocket-rescore on fpocket predictions
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
# Base directory for fpocket predictions
FPOCKET_BASE_DIR="$DATA_DIR/fpocket_preds"
# Number of parallel processes (adjust based on your system)
NUM_PROCESSES=4

# Check if p2rank executable exists
if [ ! -f "$P2RANK_BIN" ]; then
    echo "Error: p2rank executable not found at $P2RANK_BIN"
    echo "Please run get_p2rank.sh first to download and extract p2rank"
    exit 1
fi

# Counter for tracking progress
total=0
processed=0

# Iterate over each protein directory in fpocket_preds
for protein_dir in "$FPOCKET_BASE_DIR"/*; do
    if [ ! -d "$protein_dir" ]; then
        continue
    fi

    protein_id=$(basename "$protein_dir")
    
    # Find all frame directories (e.g., frame_0, frame_1, ...)
    frame_dirs=("$protein_dir"/frame_*/)
    
    # If no matching frame directories exist, skip this protein
    if [ ! -d "${frame_dirs[0]}" ]; then
        echo "No frame directories found for $protein_id, skipping."
        continue
    fi

    total=$((total + 1))
    
    dataset_file="$protein_dir/dataset.ds"

    echo "Processing protein: $protein_id"
    echo "  Writing dataset file with ${#frame_dirs[@]} frames..."

    # Create a dataset file with headers for fpocket rescoring
    {
        echo "PARAM.PREDICTION_METHOD=fpocket"
        echo ""
        echo "HEADER: prediction protein"
        echo ""
    } > "$dataset_file"

    # Append file pairs to the dataset file
    for frame_dir in "${frame_dirs[@]}"; do
        if [ -d "$frame_dir" ]; then
            # fpocket prediction pdb (e.g. frame_00000_out/frame_00000_out.pdb)
            prediction_pdb=$(find "$frame_dir" -name "*_out.pdb" | head -n 1)
            
            # Original structure pdb (e.g. frame_00000.pdb)
            frame_name=$(basename "$frame_dir" | sed 's/_out$//')
            original_pdb="$REPO_DIR/data/bioemu_results/$protein_id/${frame_name}.pdb"

            if [ -f "$prediction_pdb" ] && [ -f "$original_pdb" ]; then
                echo "$prediction_pdb $original_pdb" >> "$dataset_file"
            else
                if [ ! -f "$prediction_pdb" ]; then
                    echo "Warning: Prediction PDB not found in $frame_dir"
                fi
                if [ ! -f "$original_pdb" ]; then
                    echo "Warning: Original PDB not found at $original_pdb"
                fi
            fi
        fi
    done

    echo "  Running p2rank rescore..."
    # The output will be placed in the protein_dir by default
    "$P2RANK_BIN" rescore -o "$protein_dir" -t "$NUM_PROCESSES" "$dataset_file"
    
    if [ $? -eq 0 ]; then
        ((processed++))
        echo "  Completed"
    else
        echo "  Failed"
    fi

done

echo ""
echo "========================================="
echo "P2Rank rescore complete!"
echo "Processed: $processed/$total proteins"
echo "Results saved to respective protein folders in: $FPOCKET_BASE_DIR"
echo "========================================="
