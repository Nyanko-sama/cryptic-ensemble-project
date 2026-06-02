import os
import glob
import json
import numpy as np
import pandas as pd

def get_protein_dirs(base_dir : str, recursive=False):
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"Base directory not found: {base_dir}")

    if recursive:
        return sorted(
            [os.path.join(base_dir, name) for name in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, name))]
        )
    elif os.path.isdir(base_dir):
        return [base_dir]
    
    raise FileNotFoundError(
        f"Protein directory not found: {base_dir}."
    )

def load_all_predictions(predictions_json):
    if not os.path.isfile(predictions_json):
        raise FileNotFoundError(f"Predictions JSON file not found: {predictions_json}")
    with open(predictions_json) as f:
        all_preds = json.load(f)
    return pd.concatenate([pd.DataFrame.from_records(preds) for preds in all_preds.values()], index=all_preds.keys())

def load_evaluation_dataset(eval_dataset_path):
    df = pd.read_csv(eval_dataset_path, names=['pdb_id', 'chain_id', 'ligands', 'residue_ids'], sep=';', skipinitialspace=True)
    df['ligands'] = df['ligands'].apply(lambda x: x.split(' ') if pd.notna(x) else [])
    df['residue_ids'] = df['residue_ids'].apply(lambda x: x.split(' ') if pd.notna(x) else [])
    df['chain_id'] = df['chain_id'].str.strip()
    return df

def get_alignment_matrices(protein_dir, input_path=None):
    alignment_path = input_path or os.path.join(protein_dir, "structure_alignment_matrices.json")
    if not os.path.isfile(alignment_path):
        raise FileNotFoundError(f"Alignment matrices not found at expected path: {alignment_path}")
    with open(alignment_path) as f:
        alignments = json.load(f)["alignments"]

    for entry in alignments:
        entry["rotation"] = np.asarray(entry["rotation"], dtype=float)
        entry["translation"] = np.asarray(entry["translation"], dtype=float)
    return alignments

def gather_frames(protein_dir):
    frame_files = sorted(glob.glob(os.path.join(protein_dir, "frame_*.pdb")))
    if not frame_files:
        raise FileNotFoundError(f"No frame_*.pdb files found in protein folder: {protein_dir}")
    return frame_files

def gather_predictions(protein_dir):
    pred_files = sorted(glob.glob(os.path.join(protein_dir, "frame_*_predictions.csv")))
    if not pred_files:
        raise FileNotFoundError(f"No frame_*_predictions.csv files found in protein folder: {protein_dir}")
    
    res = []
    n_frames = len(pred_files)
    for pred_file in pred_files:
        df = pd.read_csv(pred_file, skipinitialspace=True)
        df.columns = df.columns.str.strip()  # Strip whitespace from column names
        if not {"name", "score", "center_x", "center_y", "center_z", "residue_ids"}.issubset(df.columns):
            raise ValueError(f"Missing required columns in {pred_file}. Required: name, score, center_x, center_y, center_z, residue_ids")
        
        df['frame_file'] = os.path.basename(pred_file).split('_predictions.csv')[0]
        res.append(df)

    return pd.concat(res, ignore_index=True).sort_values(by='score', ascending=False), n_frames

def gather_residues(protein_dir):
    res_files = sorted(glob.glob(os.path.join(protein_dir, "frame_*_residues.csv")))
    if not res_files:
        raise FileNotFoundError(f"No frame_*_residues.csv files found in protein folder: {protein_dir}")
    res = []
    for res in res_files:
        df = pd.read_csv(res, skipinitialspace=True)
        if not {"chain", "residue_label", "residue_name"}.issubset(df.columns):
            raise ValueError(f"Missing required columns in {res}. Required: chain, residue_label, residue_name")

        df['frame_file'] = os.path.basename(res).split('_residues.csv')[0]
        res.append(df)

    res_df = pd.concat(res, ignore_index=True).groupby(['chain', 'residue_label', 'residue_name'], sort=False).agg(list)
    return res_df