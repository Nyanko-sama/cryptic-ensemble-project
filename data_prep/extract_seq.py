'''
1. Read dataset from CryptoBench JSON file.
2. Extract pdb ids and chain ids for APO structures.
3. For each APO structure, extract sequence of specific chain and auth/label indices.
4. Save extracted sequences and indices to a new JSON file.
'''

import sys
import os
import argparse
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.working_with_cryptobench  import get_sequence_with_auth_and_label_indices, load_dataset, extract_apo_chains

parser = argparse.ArgumentParser(description='Extract sequences from cif files of APO structures and save them as JSON, with the auth ids and labels ids')
parser.add_argument('--crypto_path', type=str, default='../data/cryptobench/', help='Path to the cryptobench repository')
args = parser.parse_args()


def main(crypto_path: str):
    dataset = load_dataset(f"{crypto_path}/cryptobench-dataset/dataset.json")
    apo_chains = extract_apo_chains(dataset)

    sequences_data = {}
    for apo_pdb_id, chain_id in apo_chains.items():
        # skip proteins with multichain pocket
        if '-' in chain_id:
            continue

        sequence, label_indices, auth_indices = get_sequence_with_auth_and_label_indices(crypto_path, apo_pdb_id, chain_id)
        sequences_data[apo_pdb_id] = {
            'sequence': sequence,
            'label_indices': label_indices,
            'auth_indices': auth_indices
        }

    with open(f"{crypto_path}/cryptobench-dataset/sequences.json", 'w') as f:
        json.dump(sequences_data, f)


if __name__ == "__main__":
    main(args.crypto_path)