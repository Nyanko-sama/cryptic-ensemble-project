#!/usr/bin/env bash
set -euo pipefail

# Editable experiment configuration.
# Define each hyperparameter name and its list of values here.
# The script will run all combinations.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_SCRIPT="${PROJECT_ROOT}/src/clustering_predictor.py"
CONFORM_DIR="/auto/budejovice1/niederlj/DeepLife/bioemu_outputs/"
PREDS_DIR="${PROJECT_ROOT}/p2rank_preds"
OUTPUT_ROOT="${PROJECT_ROOT}/output"

# hyperparam_name : [list of values]
declare -A HYPERPARAMS
HYPERPARAMS[weight_type]="cosine linear none"
HYPERPARAMS[eps]="2.0 4.0 8.0"
HYPERPARAMS[score_base]="mean max median"
HYPERPARAMS[probability_agg]="mean"

# Optional fixed args for clustering_predictor.py
FIXED_ARGS=("--conform_dir" "${CONFORM_DIR}" "--preds_dir" "${PREDS_DIR}")

# If you want one run per combination only, leave DRY_RUN=false.
DRY_RUN=false

if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=true
fi

# Generate a stable ordered key list.
HYPERPARAM_KEYS=(weight_type eps score_base probability_agg)

function join_by() {
  local sep="$1"
  shift
  local out=""
  for arg in "$@"; do
    if [[ -z "${out}" ]]; then
      out="${arg}"
    else
      out+="${sep}${arg}"
    fi
  done
  printf '%s' "${out}"
}

# Recursive combination builder.
function build_combinations() {
  local idx="$1"
  shift
  local prefix=("$@")

  if [[ "${idx}" -ge "${#HYPERPARAM_KEYS[@]}" ]]; then
    echo "${prefix[*]}"
    return
  fi

  local key="${HYPERPARAM_KEYS[idx]}"
  local values=(${HYPERPARAMS[${key}]})

  for value in "${values[@]}"; do
    build_combinations "$((idx + 1))" "${prefix[@]}" "${key}=${value}"
  done
}

function run_command() {
  local cmd=("python" "${TRAIN_SCRIPT}" "${FIXED_ARGS[@]}")
  local run_name="run"

  while (( "$#" )); do
    local kv="$1"
    shift
    local key="${kv%%=*}"
    local value="${kv#*=}"
    cmd+=("--${key}" "${value}")
    run_name+="_${key}=${value}"
  done

  local output_dir="${OUTPUT_ROOT}/${run_name}"
  cmd+=("--output_dir" "${output_dir}")

  echo "========================================"
  echo "Run name: ${run_name}"
  echo "Command: ${cmd[*]}"
  echo "Output dir: ${output_dir}"
  echo "========================================"

  if [[ "${DRY_RUN}" == "true" ]]; then
    return
  fi

  mkdir -p "${output_dir}"
  printf '%s\n' "${cmd[*]}" > "${output_dir}/command.txt"
  pushd "${PROJECT_ROOT}/src" >/dev/null
  "${cmd[@]}"
  popd >/dev/null
}

# Main execution.
while IFS= read -r combo; do
  [[ -z "${combo}" ]] && continue
  IFS=' ' read -r -a entries <<< "${combo}"
  run_command "${entries[@]}"
done < <(build_combinations 0)
