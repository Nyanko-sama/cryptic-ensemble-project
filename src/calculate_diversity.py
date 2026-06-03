import argparse
import os
import sys
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
from itertools import combinations

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.data_loading import gather_frames, ordered_ca_atoms

def compute_rmsd(struct1_atoms, struct2_atoms):
    """
    Computes the RMSD between two sets of atoms.
    """
    # from Bio.PDB import Superimposer
    
    # Ensure we have the same number of atoms
    min_atoms = min(len(struct1_atoms), len(struct2_atoms))
    if min_atoms < 3:
        return np.nan # Not enough atoms to align

    super_imposer = Superimposer()
    super_imposer.set_atoms(struct1_atoms[:min_atoms], struct2_atoms[:min_atoms])
    super_imposer.apply(struct2_atoms)
    
    return super_imposer.rms

def calculate_pairwise_rmsd(protein_id, conform_dir):
    """
    Calculates the pairwise RMSD matrix for all structures of a given protein.
    """
    protein_dir = os.path.join(conform_dir, protein_id)
    frame_files = gather_frames(protein_dir)
    n_frames = len(frame_files)
    rmsd_matrix = np.zeros((n_frames, n_frames))
    
    parser = PDBParser(QUIET=True)
    
    structures = [ordered_ca_atoms(parser.get_structure("s", f)) for f in frame_files]
    
    for i in range(n_frames):
        for j in range(i + 1, n_frames):
            rmsd = compute_rmsd(structures[i], structures[j])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
            
    return pd.DataFrame(rmsd_matrix, index=[os.path.basename(f) for f in frame_files], columns=[os.path.basename(f) for f in frame_files])

def main():
    parser = argparse.ArgumentParser(description="Calculate pairwise RMSD for protein conformations.")
    parser.add_argument("--conform_dir", default='../data/bioemu_outputs', help="Directory containing conformational ensemble data.")
    parser.add_argument("--output_dir", default='../data/diversity', help="Directory to save the RMSD matrix.")
    args = parser.parse_args()

    protein_ids = [d for d in os.listdir(args.conform_dir) if os.path.isdir(os.path.join(args.conform_dir, d))]

    os.makedirs(args.output_dir, exist_ok=True)

    for protein_id in protein_ids:
        output_csv = os.path.join(args.output_dir, f"{protein_id}.csv")
        if os.path.exists(output_csv):
            print(f"Skipping {protein_id}, file already exists.")
            continue
        print(f"Processing {protein_id}...")
        
        rmsd_df = calculate_pairwise_rmsd(protein_id, args.conform_dir)
        
        rmsd_df.to_csv(output_csv)
        print(f"RMSD matrix for {protein_id} saved to {output_csv}")

if __name__ == "__main__":
    main()
