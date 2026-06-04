import argparse
import os
import sys
import json
import sys
import numpy as np
import pandas as pd

sys.path.append(os.path.join('..', "utils",))
sys.path.append(os.path.join('..', ''))

from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from utils.data_loading import gather_predictions, get_protein_dirs, gather_frames, get_alignment_matrices

def create_base_parser():   
    parser = argparse.ArgumentParser(description="Run the full prediction pipeline for cryptic pocket detection.")
    parser.add_argument("--conform_dir", default='/auto/budejovice1/niederlj/DeepLife/bioemu_outputs', help="Directory containing conformational ensemble data. Default: ../data/bioemu_results")
    parser.add_argument("--preds_dir", default="../p2rank_preds", help="Base directory where protein prediction folders are located. Defaults to ../data/p2rank_preds relative to script.")
    parser.add_argument("--output_dir", default=None, help="Directory to save the final predictions CSV file. Default: ../output")
    parser.add_argument("--ref_structure_folder", default="../data/cryptobench/cryptobench-dataset/auxiliary-data/cif-files", help="Folder containing reference structures for alignment.")
    parser.add_argument("--recursive", default=True, type=bool, help="Recursively search for protein directories under the base directory.")
    parser.add_argument("-v", type=int, default=0, help="Verbosity level. Higher values will print more detailed processing information. Default: 0 (no verbose output).")
    parser.add_argument("--alignment_dir", default="../data/prot_alignments", help="!!!REGENERATE IF YOURE USING A DIFFERENT CONFORMATION GENEARATOR!!! Path to save alignment matrices JSON file. Default is to save in each protein directory.")
    parser.add_argument("--chain_lookup_file", default="../data/chain_lookup.json", help="JSON file containing mapping of protein names to chain ids. Default: ../data/chain_lookup.json ")
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

def save_alignment_matrices(protein_dir, alignments, alignment_dir=None):
    # Extract protein name from directory path and save alignment matrices to JSON file in the same directory or specified output folder
    protein_name = os.path.basename(protein_dir)
    print(protein_name)
    output_path = os.path.join(alignment_dir, protein_name, "structure_alignment_matrices.json") if alignment_dir else os.path.join(protein_dir, "structure_alignment_matrices.json")
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
    grouped = prediction_df.groupby('frame_file')
    transform_dict = {entry['frame_file']: entry for entry in alignment_out}

    for frame, group in grouped:
        if frame not in transform_dict:
            raise ValueError(f"No alignment information found for frame: {frame}")
        R = transform_dict[frame]['rotation']
        t = transform_dict[frame]['translation']
        coords = group[["center_x", "center_y", "center_z"]].values
        aligned_coords = np.dot(np.vstack(coords), R.T) + t
        prediction_df.loc[group.index, ["center_x", "center_y", "center_z"]] = aligned_coords

    return prediction_df

def compute_structure_alignments(protein_dir : str, reference_structure_path : str, verbose=False, alignment_dir=None):
    # Check if alignment matrices already exist
    prot_name = os.path.basename(protein_dir)
    final_align_dir = os.path.join(alignment_dir, prot_name) if alignment_dir else protein_dir
    if os.path.isfile(os.path.join(final_align_dir, "structure_alignment_matrices.json")):
        if verbose:
            print(f"Alignment matrices already exist for {prot_name}. Loading from file.")
        return get_alignment_matrices(final_align_dir)

    frame_files = gather_frames(protein_dir)
    if verbose:
        print(f"Found {len(frame_files)} frame files for {protein_dir}. Computing alignments to reference structure at: {reference_structure_path}")
    
    # There was a bug MMCIF parser that caused it to fail when trying to warn about missing residues. Had to manually fix the 
    # parameters inside the source code of the parser to get it working. See Bio/PDB/MMCIFParser.py line 241 if something like that happens again.
    parser = MMCIFParser(QUIET=True, auth_residues=False)

    original_structure = parser.get_structure(os.path.basename(reference_structure_path), reference_structure_path)
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

    save_alignment_matrices(protein_dir, alignments, alignment_dir=alignment_dir)
    return get_alignment_matrices(final_align_dir)

def load_chain_lookup(path):
    """Load chain information from a JSON file."""
    with open(path, 'r') as f:
        chain_lookup = json.load(f)
    return chain_lookup

def prediction_pipeline(args, aggregation_func, process_func, output_path):
    # Load chain information from both train and test CSV files
    chain_lookup = load_chain_lookup(args.chain_lookup_file)
    all_preds = {}
    
    for protein_dir in get_protein_dirs(args.preds_dir, recursive=args.recursive):
        prot_name = os.path.basename(protein_dir)

        if args.v > 0:
            print(f"Processing protein: {prot_name}")

        try:
            predictions_df, n_frames = gather_predictions(protein_dir)
        except Exception as e:
            print(f"Error occurred while gathering predictions for {prot_name}: {e}")

        if predictions_df.empty:
            print(f"No predictions in the prediction files for {prot_name}. Skipping.")
            continue
        reference_structure_path = os.path.join(args.ref_structure_folder, prot_name + ".cif")
        if args.v > 0:
            print(f"Computing structure alignments for {prot_name} using reference structure at: {reference_structure_path}")

        # Align pocket coordinates to reference structure frame
        alignment_out = compute_structure_alignments(os.path.join(args.conform_dir, prot_name), reference_structure_path=reference_structure_path, verbose=args.v > 0,
                                                     alignment_dir=args.alignment_dir)
        predictions_df = align_pocket_coordinates(predictions_df, alignment_out)

        # Aggregate pockets and process final predictions for this protein
        aggregation_out = aggregation_func(predictions_df, verbose=args.v)
        final_pred_df = process_func(aggregation_out, n_frames=n_frames, verbose=args.v)

        # Add chain column from lookup
        if prot_name in chain_lookup:
            final_pred_df['chain'] = chain_lookup[prot_name]
        else:
            raise ValueError(f"Chain information not found for {prot_name} in CSV files.")

        if args.v > 0:
            print(f"Finished processing {prot_name}. Saving final predictions to: {output_path}")

        # Save final predictions for this protein and add to overall results
        all_preds[os.path.basename(protein_dir)] = final_pred_df.to_dict(orient='records')
        save_predictions(final_pred_df, os.path.join(output_path, os.path.basename(protein_dir) + "_aggregated_predictions.csv"))

    save_all_predictions(all_preds, output_path, description=os.path.basename(output_path))

def save_all_predictions(all_preds, output_path, team_name="Praga", model_version="v1.0", 
                         submission_date=None, description="Cryptic pocket predictions using ensemble methods"):
    """Save predictions in the standardized JSON format.
    
    Args:
        all_preds: Dictionary mapping protein names to lists of prediction records
        output_path: Path to save the JSON file
        team_name: Team name for metadata
        model_version: Model version for metadata
        submission_date: Submission date (default: today's date)
        description: Description of the submission
    """
    from datetime import datetime
    
    if submission_date is None:
        submission_date = datetime.now().strftime("%Y-%m-%d")
    
    # Build predictions array
    predictions = []
    for pdb_id, pred_records in all_preds.items():
        if not pred_records:
            continue
        # Get chain and protein-level info from first record
        first_record = pred_records[0]
        chain = first_record.get('chain')
        
        # Build ranked pockets by sorting by score/probability
        ranked_pockets = []
        for rank, record in enumerate(sorted(pred_records, 
                                            key=lambda x: x.get('score'), 
                                            reverse=True), 1):
            pocket = {
                "rank": rank,
                "probability" : float(record.get('probability')),
                "score": float(record.get('score')), # For now, we are using the score as the probability. This can be changed if a separate probability field is available.
                "residues":[f"{chain}:{residue[2:]}" for residue in record.get('residue_ids', [])],
                "center": [
                    float(record.get('center_x', 0)),
                    float(record.get('center_y', 0)),
                    float(record.get('center_z', 0))
                ]
            }
            ranked_pockets.append(pocket)
        
        predictions.append({
            "pdb_id": pdb_id,
            "chain": chain,
            "ranked_pockets": ranked_pockets
        })
    
    # Build final JSON structure
    output_json = {
        "metadata": {
            "team_name": team_name,
            "model_version": model_version,
            "submission_date": submission_date,
            "description": description
        },
        "predictions": predictions
    }
    
    # Save to file
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    
    with open(os.path.join(output_path, "all_predictions.json"), 'w') as f:
        json.dump(output_json, f, indent=2)

def save_predictions(prediction_df, output_path):
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    prediction_df.to_csv(output_path, index=False)

def baseline_eval_file_create(args, output_path):
    chain_lookup = load_chain_lookup(args.chain_lookup_file)
    all_preds = {}

    for protein_dir in get_protein_dirs(args.preds_dir, recursive=args.recursive):
        prot_name = os.path.basename(protein_dir)

        if args.v > 0:
            print(f"Processing protein: {prot_name}")

        try:
            predictions_df, n_frames = gather_predictions(protein_dir)
        except Exception as e:
            print(f"Error occurred while gathering predictions for {prot_name}: {e}")
            continue

        if predictions_df.empty:
            print(f"No predictions in the prediction files for {prot_name}. Skipping.")
            continue
        
        # Add chain column from lookup
        if prot_name in chain_lookup:
            predictions_df['chain'] = chain_lookup[prot_name]
        else:
            raise ValueError(f"Chain information not found for {prot_name} in CSV files.")
        
        all_preds[prot_name] = predictions_df.to_dict(orient='records')
        
        if args.v > 0:
            print(f"Finished processing {prot_name}. Saving final predictions to: {output_path}")
     
        valid_records = []
        for record in all_preds[prot_name]:
            residue_ids = record.get('residue_ids')
            if isinstance(residue_ids, str):
                record['residue_ids'] = residue_ids.split()
                valid_records.append(record)
            else:
                if args.v > 0:
                    print(f"Skipping record with unexpected format for residue_ids: {record}")
        all_preds[prot_name] = valid_records

    save_all_predictions(all_preds, output_path, description=os.path.basename(output_path))