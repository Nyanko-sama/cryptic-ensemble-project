import os
import glob
import pandas as pd

def get_protein_dirs(protein_arg : str, base_dir : str, recursive=False):
    candidate = os.path.join(base_dir, protein_arg if protein_arg else "")
    if os.path.isdir(candidate):
        if recursive:
            return sorted(
                [os.path.join(candidate, name) for name in os.listdir(candidate) if os.path.isdir(os.path.join(candidate, name))]
            )
        return [candidate]

    if os.path.isdir(protein_arg):
        return [protein_arg]

    raise FileNotFoundError(
        f"Protein directory not found: {protein_arg}. Expected under base dir: {base_dir}."
    )

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
        
        df['frame'] = os.path.basename(pred_file).split('_predictions.csv')[0]
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

        df['frame'] = os.path.basename(res).split('_residues.csv')[0]
        res.append(df)

    res_df = pd.concat(res, ignore_index=True).groupby(['chain', 'residue_label', 'residue_name'], sort=False).agg(list)
    return res_df