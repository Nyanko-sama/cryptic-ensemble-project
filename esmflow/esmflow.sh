#!/bin/bash
#PBS -N esmflow
#PBS -l select=1:ncpus=1:ngpus=1:gpu_mem=24gb:mem=32gb:scratch_local=20gb
#PBS -l walltime=48:00:00
#PBS -M neliakorotushak@gmail.com
#PBS -o esmflow.txt
#PBS -j oe

# ==========================================
# 1. SCRATCHDIR & SINGULARITY SETUP
# ==========================================

# Fail-safe: Check if scratch directory was successfully allocated
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# Force Singularity to use the compute node's scratch space for all temp/cache files.
# This prevents you from running out of space in your Metacentrum home directory.
export SINGULARITY_CACHEDIR=$SCRATCHDIR
export SINGULARITY_LOCALCACHEDIR=$SCRATCHDIR
export SINGULARITY_TMPDIR=$SCRATCHDIR

# ==========================================
# 2. RUNNING THE CONTAINER
# ==========================================

# NOTE: Update the image tag below (currently esmflow-container:v1) to whatever 
# name/tag you used when pushing your new Docker image to Docker Hub!
singularity exec --nv --bind /storage docker://neliakorotushak/esmflow:v1 bash -c '

# ----------------------------------------------------
# EVERYTHING BELOW THIS LINE RUNS INSIDE THE CONTAINER
# ----------------------------------------------------

# Internal container environment setup
CONTAINER_PYTHON="/opt/mambaforge/envs/alphaflow/bin/python"
CONTAINER_ALPHAFLOW_DIR="/opt/workspace/alphaflow"
export PYTHONPATH="${CONTAINER_ALPHAFLOW_DIR}:${PYTHONPATH}"

# Your exact Metacentrum directories
MY_DIR="/storage/praha1/home/nelia_k/alphaflow"
INPUT_CSV="$MY_DIR/input.csv"
FINAL_RESULTS_DIR="$MY_DIR/results"

echo "=== STEP 1: Running ESMFlow Inference ==="

# Note: ESMFold is a single-sequence model, so it bypasses MSAs. 
$CONTAINER_PYTHON $CONTAINER_ALPHAFLOW_DIR/predict.py \
    --mode esmfold \
    --input_csv $INPUT_CSV \
    --weights /opt/workspace/alphaflow/params/esmflow_pdb_base_202402.pt \
    --samples 40 \
    --outpdb $FINAL_RESULTS_DIR

echo "=== Pipeline Completed Successfully ==="
'