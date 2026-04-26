'''
1. Read sequences from JSON file. 
2. For each sequence, run BioEMU to generate samples.
3. Convert BioEMU trajectory to PDB files and set residue IDs based on auth_indices from the JSON.
4. Save PDB files to output directory.
'''


import argparse
import numpy as np
from bioemu.sample import main as sample
import os
import csv 
import MDAnalysis as mda
import pandas as pd
import json

parser = argparse.ArgumentParser(description='Run BioEMU convergence analysis.')

parser.add_argument('--n', type=int, default=5, help='Number of BioEMU samples to test')
parser.add_argument('--sequences', type=str, default="/home/nelia/deeplife/data_prep/20_train_sequences.json", help='Path to JSON file containing sequences')
parser.add_argument('--output_dir', type=str, default="/home/nelia/deeplife/data/bioemu", help='Directory to save results')


def bioemu_to_pdb(bioemu_path: str, residue_indices: list[int]) -> list[str]:
    """Convert BioEMU trajectory to PDB files and set residue IDs."""

    topology_path = os.path.join(bioemu_path, "topology.pdb")
    xtc_path = os.path.join(bioemu_path, "samples.xtc")

    u = mda.Universe(topology_path, xtc_path)
    protein = u.select_atoms("all")

    # set residue indices
    if len(protein.residues) != len(residue_indices):
        raise ValueError("Length of residue_indices does not match the number of residues.")
    protein.residues.resids = np.array(residue_indices)

    proteins_paths = []

    for i, ts in enumerate(u.trajectory):
        out_path = os.path.join(bioemu_path, f"frame_{i:05d}.pdb")
        protein.write(out_path)
        proteins_paths.append(out_path)

    return proteins_paths



def run_bioemu_all_sequences(n:int, sequences_list:dict, output_dir: str) -> bool:

    os.makedirs(f"{output_dir}", exist_ok=True)

    for protein_id, data in sequences_list.items():
        path_bioemu = f"{output_dir}/{protein_id}"
        os.makedirs(path_bioemu, exist_ok=True)

        seq = data['sequence']

        # Retry once on known BioEMU suffix error; skip protein if it persists.
        success = False
        for attempt in range(2):
            try:
                sample(sequence=seq, num_samples=n, output_dir=path_bioemu)
                success = True
                break
            except ValueError as exc:
                if "Invalid suffix" in str(exc) and "_unphysical.xtc" in str(exc):
                    if attempt == 0:
                        print(f"Retrying BioEMU for {protein_id} due to unphysical samples.")
                        continue
                    print(f"Skipping {protein_id} after repeated unphysical samples error.")
                    success = False
                else:
                    raise

        if not success:
            continue
        
        # get auth_indices list and convert all elements to int
        residue_indices = [int(x) for x in data['auth_indices']]
        path_pdbs = bioemu_to_pdb(path_bioemu, residue_indices)
        if len(path_pdbs) < 2:
            print(f"Skipping {protein_id}: fewer than 2 valid proteins after filtering.")
            continue

        break

def main():
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.sequences) as f:
        sequences = json.load(f)

    run_bioemu_all_sequences(args.n, sequences, args.output_dir)


if __name__ == "__main__":
    main()



