#!/bin/bash
#PBS -N deeplife_bioemu
#PBS -l select=1:ncpus=8:mem=64gb:ngpus=1:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -j oe

set -euo pipefail

START_TIME=$(date +%s)

############################################
# USER CONFIG
############################################
# Folder with your python files to copy to scratch (must include your main script + rmsd_convergence.py)
PROJECT_DIR="/storage/praha1/home/nelia_k/deeplife"

# Your main script filename inside PROJECT_DIR
MAIN_PY="run_bioemu.py"

# Output dir (final persistent storage)
OUTDIR="${PROJECT_DIR}/bioemu_results"

# Optional: sequences.json path used by your script
SEQUENCES_JSON="${PROJECT_DIR}/20_train_sequences.json"

# Conda env name (stored in your home by default)
ENV_NAME="bioemu_gpu"


############################################
# MODULES + GPU CHECK
############################################
module purge
module load mambaforge

echo "=== Host: $(hostname) ==="
echo "=== PBS jobid: ${PBS_JOBID:-unknown} ==="

# Hard-stop if GPU is not present/usable
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "ERROR: nvidia-smi not found -> no GPU environment on this node. Exiting."
  exit 1
fi

if ! nvidia-smi -L >/dev/null 2>&1; then
  echo "ERROR: GPU not accessible (nvidia-smi failed). Exiting."
  nvidia-smi || true
  exit 1
fi

echo "=== GPUs detected ==="
nvidia-smi -L

############################################

############################################
# SCRATCH SETUP (MetaCentrum native)
############################################

if [[ -z "${SCRATCHDIR:-}" ]]; then
  echo "ERROR: SCRATCHDIR not set. Did you request scratch_local?"
  exit 1
fi

echo "=== Using SCRATCHDIR: $SCRATCHDIR ==="
echo "=== Scratch type: ${SCRATCH_TYPE:-unknown} ==="
echo "=== Scratch volume: ${SCRATCH_VOLUME:-unknown} ==="

WORKDIR="$SCRATCHDIR/work"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Copy your scripts to scratch
echo "=== Copying project to scratch ==="
rsync -av --delete "${PROJECT_DIR}/" "$WORKDIR/"

# Copy sequences.json to scratch if it exists
if [[ -f "$SEQUENCES_JSON" ]]; then
  cp "$SEQUENCES_JSON" "$WORKDIR/sequences.json"
  SEQUENCES_JSON="$WORKDIR/sequences.json"
  echo "Copied sequences.json to scratch: $SEQUENCES_JSON"
else
  echo "WARNING: SEQUENCES_JSON not found at $SEQUENCES_JSON"
fi

# Ensure output dir exists (persistent)
mkdir -p "$OUTDIR"


############################################
# ACTIVATE EXISTING ENV
############################################
ENV_PREFIX="/storage/praha1/home/nelia_k/conda_envs/bioemu_gpu"

echo "=== Activating environment at $ENV_PREFIX ==="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_PREFIX"

# Verify environment is activated
echo "Active environment: $CONDA_DEFAULT_ENV"
echo "Environment python: $(which python)"

# caches on scratch (safe)
export HF_HOME="$SCRATCHDIR/hf_home"
export TRANSFORMERS_CACHE="$SCRATCHDIR/hf_cache"
export XDG_CACHE_HOME="$SCRATCHDIR/xdg_cache"
export TMPDIR="$SCRATCHDIR/tmp"
mkdir -p "$HF_HOME" "$TRANSFORMERS_CACHE" "$XDG_CACHE_HOME" "$TMPDIR"

echo "=== Python ==="
python -V
echo "=== Pip ==="
pip -V

############################################
# Sanity checks
############################################
echo "=== Torch CUDA check ==="
python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("gpu:", torch.cuda.get_device_name(0))
PY

############################################

############################################
# RUN YOUR SCRIPT
############################################
echo "=== Running convergence analysis ==="
echo "Project files in: $WORKDIR"
echo "Output dir: $OUTDIR"
echo "Sequences file: $SEQUENCES_JSON"

# Run with args (adjust n_min/n_max as needed)
python "$WORKDIR/$MAIN_PY" \
  --n 20 \
  --output_dir "$OUTDIR" \
  --sequences "$SEQUENCES_JSON"

echo "=== Script finished ==="

############################################
# SAVE LOGS / ARTIFACTS BACK
############################################
# Save job log & a snapshot of scratch work (optional)
JOBLOG_NAME="job_${PBS_JOBID:-$$}.log"
echo "Saving scratch snapshot + logs..."

# Save a copy of stdout/stderr (PBS -j oe usually puts it in the submission dir;
# we also capture a run log in OUTDIR)
{
  echo "=== Job metadata ==="
  date
  echo "Host: $(hostname)"
  echo "JobID: ${PBS_JOBID:-unknown}"
  echo "SCRATCHDIR: $SCRATCHDIR"
  echo
  echo "=== GPU ==="
  nvidia-smi || true
} > "${OUTDIR}/${JOBLOG_NAME}"

# Optional: copy any additional files produced in scratch project folder
# (Your script writes directly to OUTDIR, so this is mostly for debugging)
rsync -av --exclude="__pycache__" "$WORKDIR/" "${OUTDIR}/scratch_debug_${PBS_JOBID:-$$}/" || true

############################################
# CLEANUP
############################################
echo "=== Cleaning up scratch ==="
cd /
rm -rf "$SCRATCHDIR"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Script finished in $(date -d@$RUNTIME -u +%H:%M:%S)"

echo "=== DONE ==="