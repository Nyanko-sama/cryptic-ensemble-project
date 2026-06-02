import biotite.structure.io.pdbx as pdbx
from biotite.structure.io.pdbx import get_structure
from biotite.structure import get_residues
import numpy as np
import pandas as pd
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', ''))

from utils.data_loading import load_cryptobench_csv


def compute_center(points):
    return np.mean(points, axis=0)

def extract_site_coords(cif_file_path, chain_ids, binding_residues, auth=True):
    cif_file = pdbx.CIFFile.read(cif_file_path)
    
    protein = get_structure(cif_file, model=1, use_author_fields=auth)
    res = []
    for chain_id in chain_ids:   
        protein = protein[(protein.atom_name == "CA") 
                            & (protein.element == "C") 
                            & (protein.chain_id == chain_id) ]

        residue_ids, residue_types = get_residues(protein)
        coordinates = []
        for i in range(len(residue_ids)):
            if residue_ids[i] in binding_residues:
                coordinates.append(protein[i].coord)

        res.append(coordinates)
    return res

def extract_ids(cif_file_path, chain_id, auth=True):
    import biotite.structure.io.pdbx as pdbx
    from biotite.structure.io.pdbx import get_structure
    from biotite.structure import get_residues

    cif_file = pdbx.CIFFile.read(cif_file_path)
    
    protein = get_structure(cif_file, model=1, use_author_fields=auth)
    protein = protein[(protein.atom_name == "CA") 
                        & (protein.element == "C") 
                        & (protein.chain_id == chain_id) ]

    residue_ids, residue_types = get_residues(protein)
    return residue_ids

def create_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Precompute evaluation data for Cryptobench")
    parser.add_argument("--csv_path", type=str, required=True, help="Path to a CSV file containing columns: pdb_id, chain_id, ligands, residue_ids")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save precomputed evaluation data (numpy file)")
    parser.add_argument("--auth_labels", action='store_true', help="Whether to use author residue numbering in labels")
    return parser

def precompute_addition_eval_data(dataframe, auth_labels=True):
    centers = []
    ids = []
    for idx, row in dataframe.iterrows():
        print(f"Processing {idx}...")
        cif_file_path = f"../data/cryptobench/cryptobench-dataset/auxiliary-data/cif-files/{idx.lower()}.cif"
        coords = [extract_site_coords(cif_file_path, row['chain_id'], binding_residues, auth=auth_labels) for binding_residues in row['residue_ids']]
        temp_ids = [extract_ids(cif_file_path, chain_id, auth=auth_labels).tolist() for chain_id in row['chain_id']]
        center = [compute_center(np.vstack(coord)).tolist() for coord in coords]
        centers.append(center)
        ids.append(temp_ids)
    return centers, ids

def prepare_dataset(args):
    eval_dataset = load_cryptobench_csv(args.csv_path).groupby('pdb_id').agg(list)
    print(eval_dataset.head(1))
    centers, ids = precompute_addition_eval_data(eval_dataset, auth_labels=args.auth_labels)
    eval_dataset['centers'] = centers
    eval_dataset['total_residue_ids'] = ids
    eval_dataset.rename(columns={'residue_ids': 'binding_residue_ids'}, inplace=True)
    output_path = f"{args.output_dir}/{os.path.basename(args.csv_path).replace('.csv', '')}_eval_dataset_{'auth' if args.auth_labels else 'non_auth'}_labels.json"
    eval_dataset.to_json(output_path, indent=2)

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    prepare_dataset(args)