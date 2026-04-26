#!/bin/bash

# Set the root directory of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
DATA_DIR="${REPO_ROOT}/data"
CRYPTOBENCH_PATH="${DATA_DIR}/cryptobench"

# Create the data directory in the root of the repository
mkdir -p "${DATA_DIR}"

# Change to the data directory
cd "${DATA_DIR}" || exit

# Download and prepare the CryptoBench dataset if it doesn't exist
if [ ! -d "cryptobench" ]; then
    wget -O cryptobench.zip https://files.de-1.osf.io/v1/resources/pz4a9/providers/osfstorage/?zip= --no-check-certificate
    unzip cryptobench.zip
    rm cryptobench.zip
fi

cd "${REPO_ROOT}/data_prep" || exit

# Unzip the cif-files.zip and remove the zip file
unzip "${CRYPTOBENCH_PATH}/cryptobench-dataset/auxiliary-data/cif-files.zip" -d "${CRYPTOBENCH_PATH}/cryptobench-dataset/auxiliary-data"
rm "${CRYPTOBENCH_PATH}/cryptobench-dataset/auxiliary-data/cif-files.zip"


# Extract sequences from the cif files of APO structure and save them as JSON, with the auth ids and labels ids
python extract_seq.py --crypto_path "${CRYPTOBENCH_PATH}"
