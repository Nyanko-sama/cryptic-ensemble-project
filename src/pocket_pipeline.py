
import argparse
import os
import sys
import json
import sys
import numpy as np

sys.path.append(os.path.join('..', "utils",))
sys.path.append(os.path.join('..', ''))

from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from utils.data_loading import gather_predictions, get_protein_dirs, gather_frames, get_alignment_matrices

def create_base_parser():   
    parser = argparse.ArgumentParser(description="Run the full prediction pipeline for cryptic pocket detection.")
    parser.add_argument("--conform_dir", default='../data/bioemu_results', help="Directory containing conformational ensemble data. Default: ../data/bioemu_results")
    parser.add_argument("--preds_dir", default="../data/p2rank_preds", help="Base directory where protein prediction folders are located. Defaults to ../data/p2rank_preds relative to script.")
    parser.add_argument("--output_path", default="final_predictions.csv", help="Path to save the final predictions CSV file. Default: final_predictions.csv")
    parser.add_argument("--ref_structure_folder", default="../data/cryptobench/cryptobench-dataset/auxiliary-data/cif-files", help="Folder containing reference structures for alignment.")
    parser.add_argument("--recursive", default=True, type=bool, help="Recursively search for protein directories under the base directory.")
    parser.add_argument("-v", type=int, default=0, help="Verbosity level. Higher values will print more detailed processing information. Default: 0 (no verbose output).")
    return parser


def _ordered_ca_atoms(structure):
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    atoms.append((chain.id, residue.get_id()[1], residue["CA"]))
        break
    atoms.sort(key=lambda item: (item[0], item[1]))
    return [entry[2] for entry in atoms]

def save_alignment_matrices(protein_dir, alignments, output_path=None):
    output_path = output_path or os.path.join(protein_dir, "structure_alignment_matrices.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump({
            "alignments": alignments,
        }, f, indent=2)
    return output_path

def compute_frame_alignment(ref_atoms, frame_atoms):
    n_matched = min(len(ref_atoms), len(frame_atoms))
    if n_matched < 3:
        raise ValueError("At least 3 Cα atoms are required for alignment.")
    sup = Superimposer()
    sup.set_atoms(ref_atoms[:n_matched], frame_atoms[:n_matched])
    R = np.asarray(sup.rotran[0], dtype=float)
    t = np.asarray(sup.rotran[1], dtype=float)
    rmsd = float(sup.rms)
    return R, t, n_matched, rmsd

def align_pocket_coordinates(prediction_df, alignment_out):
    grouped = prediction_df.groupby('frame')
    transform_dict = {entry['frame']: entry for entry in alignment_out}

    for frame, group in grouped:
        if frame not in transform_dict:
            raise ValueError(f"No alignment information found for frame: {frame}")
        R = transform_dict[frame]['rotation']
        t = transform_dict[frame]['translation']
        coords = group[["center_x", "center_y", "center_z"]].values
        aligned_coords = np.dot(np.vstack(coords), R.T) + t
        prediction_df.loc[group.index, ["center_x", "center_y", "center_z"]] = aligned_coords

    return prediction_df

def compute_structure_alignments(protein_dir : str, reference_structure_path : str, verbose=False):
    # Check if alignment matrices already exist
    if os.path.isfile(os.path.join(protein_dir, "structure_alignment_matrices.json")):
        if verbose:
            print(f"Alignment matrices already exist for {protein_dir}. Loading from file.")
        return get_alignment_matrices(protein_dir)

    frame_files = gather_frames(protein_dir)
    if verbose:
        print(f"Found {len(frame_files)} frame files for {protein_dir}. Computing alignments to reference structure at: {reference_structure_path}")
    
    # There was a bug MMCIF parser that caused it to fail when trying to warn about missing residues. Had to manually fix the 
    # parameters inside the source code of the parser to get it working. See Bio/PDB/MMCIFParser.py line 241 if something like that happens again.
    parser = MMCIFParser(QUIET=True, auth_residues=False)
    original_structure = parser.get_structure(os.path.basename(reference_structure_path).split('.')[0], reference_structure_path)
    reference_atoms = _ordered_ca_atoms(original_structure)

    alignments = []
    for frame_file in frame_files:
        frame_structure = PDBParser(QUIET=True).get_structure("frame", frame_file)
        frame_atoms = _ordered_ca_atoms(frame_structure)
        R, t, n_matched, rmsd = compute_frame_alignment(reference_atoms, frame_atoms)
        alignments.append({
            "frame_file": os.path.basename(frame_file),
            "matched_ca_atoms": int(n_matched),
            "rotation": R.tolist(),
            "translation": t.tolist(),
            "rmsd": rmsd
        })

    save_alignment_matrices(protein_dir, alignments)
    return get_alignment_matrices(protein_dir)

def prediction_pipeline(args, aggregation_func, process_func, output_path):
    all_preds = {}
    for protein_dir in get_protein_dirs(args.preds_dir, recursive=args.recursive):
        prot_name = os.path.basename(protein_dir)

        if args.v > 0:
            print(f"Processing protein: {prot_name}")

        predictions_df = gather_predictions(protein_dir)
        if predictions_df.empty:
            raise ValueError(f"No predictions found in protein directory: {protein_dir}")
        n_frames = predictions_df.shape[0]
        reference_structure_path = os.path.join(args.ref_structure_folder, prot_name + ".cif")
        if args.v > 0:
            print(f"Computing structure alignments for {prot_name} using reference structure at: {reference_structure_path}")

        # Align pocket coordinates to reference structure frame
        alignment_out = compute_structure_alignments(os.path.join(args.conform_dir, prot_name), reference_structure_path=reference_structure_path, verbose=args.v > 0)
        predictions_df = align_pocket_coordinates(predictions_df, alignment_out)

        # Aggregate pockets and process final predictions
        aggregation_out = aggregation_func(predictions_df, verbose=args.v)
        final_pred_df = process_func(aggregation_out, n_frames=n_frames, verbose=args.v)

        if args.v > 0:
            print(f"Finished processing {prot_name}. Saving final predictions to: {output_path}")

        # Save final predictions for this protein and add to overall results
        all_preds[os.path.basename(protein_dir)] = final_pred_df.to_dict(orient='records')
        save_predictions(final_pred_df, os.path.join(output_path, os.path.basename(protein_dir) + "_final_predictions.csv"))

    save_all_predictions(all_preds, output_path)
    
def save_all_predictions(all_preds, output_path):
    # Save all predictions to a JSON file for easier downstream analysis
    with open(os.path.join(output_path, "all_predictions.json"), 'w') as f:
        json.dump(all_preds, f, indent=4)

def save_predictions(prediction_df, output_path):
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    prediction_df.to_csv(output_path, index=False)