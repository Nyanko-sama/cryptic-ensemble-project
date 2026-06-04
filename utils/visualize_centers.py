#!/usr/bin/env python3
"""
visualize.py — DeepLife 2026 PyMOL Visualization Generator

Generates a PyMOL script (.pml) to visually inspect cryptic pocket predictions
for a specific target against the ground truth.

Usage:
    python visualize_centers.py \
        --predictions src/all_predictions.json \
        --ground-truth data/test.csv \
        --structures data/structures/ \
        --target 1abc \
        --output visualize_1abc.pml


python visualize_centers.py --predictions ../output/fpocket_eps_4/all_predictions.json --ground-truth ../data/test.csv --target 1arl --output ../visualization/visualize_1arl.pml --structures ../data/cryptobench/cryptobench-dataset/auxiliary-data/cif-files/
"""

import argparse
import json
import sys
from pathlib import Path
import numpy as np
import biotite.structure.io.pdbx as pdbx
from biotite.structure.io.pdbx import get_structure
from biotite.structure import get_residues

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_ground_truth(csv_path):
    """
    Returns dict: (pdb_id_lower, chain) -> list[list[str]]
    Each inner list is the auth_seq_id residues for one pocket as strings.
    """
    gt = {}
    with open(csv_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(";")
            pdb_id = parts[0].lower()
            chain = parts[1]
            residues = parts[3].split()
            gt.setdefault((pdb_id, chain), []).append(residues)
    return gt

def load_submission(json_path):
    with open(json_path) as f:
        return json.load(f)

# ---------------------------------------------------------------------------
# Residue parsing
# ---------------------------------------------------------------------------

def parse_residues(residue_list):
    """
    Strips chain prefix if present; preserves insertion codes.
    Handles both 'A:220' and '220' cleanly.
    """
    return [r.split(":", 1)[-1] for r in residue_list]

def to_pymol_selection(name, chain, residues):
    """
    Generates a PyMOL selection string. 
    Groups residues using the '+' operator.
    """
    clean_res = parse_residues(residues)
    resi_str = "+".join(clean_res)
    return f"select {name}, chain {chain} and resi {resi_str}"

# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def _load_ca_coords(cif_path, chain_id):
    """Returns dict {residue_key: np.array([x,y,z])} for all Cα in chain."""
    cif_file = pdbx.CIFFile.read(cif_path)
    protein = get_structure(cif_file, model=1, use_author_fields=True)
    protein = protein[
        (protein.atom_name == "CA")
        & (protein.element == "C")
        & (protein.chain_id == chain_id)
    ]
    residue_ids, _ = get_residues(protein)

    coords_all = protein.coord
    coords = {}
    for i, atom in enumerate(protein):
        ins = atom.ins_code.strip() if hasattr(atom, "ins_code") else ""
        key = str(atom.res_id) + ins
        coords[key] = coords_all[i]

    return coords

def _centroid(ca_coords, residue_keys):
    clean_keys = parse_residues(residue_keys)
    pts = [ca_coords[r] for r in clean_keys if r in ca_coords]
    if not pts:
        raise ValueError(f"No Cα atoms found for residues {clean_keys}")
    return np.mean(pts, axis=0)

def _get_ca_cache(pdb_id, chain, ca_coords_cache, structures_dir):
    key = (pdb_id.lower(), chain)
    if key not in ca_coords_cache:
        cif_path = Path(structures_dir) / f"{pdb_id.lower()}.cif"
        if not cif_path.exists():
            raise FileNotFoundError(f"Structure file not found: {cif_path}")
        ca_coords_cache[key] = _load_ca_coords(cif_path, chain)
    return ca_coords_cache[key]

def resolve_center(pocket, pdb_id, chain, ca_coords_cache, structures_dir):
    """Returns np.array [x,y,z]; uses pocket['center'] if present."""
    center = pocket.get("center")
    if center is not None:
        return np.array(center, dtype=float)
    ca = _get_ca_cache(pdb_id, chain, ca_coords_cache, structures_dir)
    return _centroid(ca, pocket["residues"])

def gt_center(residue_keys, pdb_id, chain, ca_coords_cache, structures_dir):
    ca = _get_ca_cache(pdb_id, chain, ca_coords_cache, structures_dir)
    return _centroid(ca, residue_keys)

# ---------------------------------------------------------------------------
# CLI & Main Generation
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Generate PyMOL script for Top-N cryptic pockets.")
    p.add_argument("--predictions", required=True, help="Path to predictions.json")
    p.add_argument("--ground-truth", required=True, help="Path to test.csv")
    p.add_argument("--structures", required=True, help="Directory with .cif files")
    p.add_argument("--target", required=True, help="Target ID as PDB or PDB:CHAIN (e.g., 1abc or 1abc:A)")
    p.add_argument("--output", default="visualize.pml", help="Output PyMOL script (.pml)")
    return p.parse_args()

def main():
    args = parse_args()
    
    # Parse target string
    if ":" in args.target:
        target_pdb, target_chain = args.target.split(":", 1)
        target_pdb = target_pdb.lower()
    else:
        target_pdb = args.target.lower()
        target_chain = None

    submission = load_submission(args.predictions)

    # Infer chain if not provided
    if target_chain is None:
        matches = [p for p in submission["predictions"] if p["pdb_id"].lower() == target_pdb]
        if not matches:
            sys.exit(f"Error: PDB '{target_pdb}' not found in {args.predictions}")
        
        target_chain = matches[0]["chain"]
        if len(matches) > 1:
            print(f"Note: Multiple chains found for '{target_pdb}'. Auto-selected chain '{target_chain}'. "
                  f"To specify another, use --target {target_pdb}:CHAIN")
        else:
            print(f"Inferred chain '{target_chain}' for PDB '{target_pdb}'.")

    # Locate prediction entry
    pred_entry = next((p for p in submission["predictions"] 
                       if p["pdb_id"].lower() == target_pdb and p["chain"] == target_chain), None)

    if not pred_entry:
        sys.exit(f"Error: Target {target_pdb}:{target_chain} not found in {args.predictions}")

    # Load ground truth
    gt = load_ground_truth(args.ground_truth)
    key = (target_pdb, target_chain)
    
    if key not in gt:
        sys.exit(f"Error: Target {target_pdb}:{target_chain} not found in ground truth ({args.ground_truth})")

    true_pockets = gt[key]
    N_true = 2
    print(f"Target {target_pdb}:{target_chain} has {N_true} true pocket(s).")

    ranked = sorted(pred_entry["ranked_pockets"], key=lambda p: p["rank"])
    top_n_preds = ranked[:N_true]

    ca_coords_cache = {}
    cif_path = Path(args.structures) / f"{target_pdb}.cif"

    pml_lines = []
    
    # a) experimental protein
    pml_lines.append(f"load {cif_path.absolute()}, protein")
    pml_lines.append("hide everything, protein")
    pml_lines.append("show cartoon, protein")
    pml_lines.append("color white, protein")
    pml_lines.append("")

    # True Pockets processing
    pml_lines.append("# --- GROUND TRUTH POCKETS ---")
    for i, tp in enumerate(true_pockets):
        center = gt_center(tp, target_pdb, target_chain, ca_coords_cache, args.structures)
        
        # b) real center as a sphere
        pml_lines.append(f"pseudoatom true_center_{i+1}, pos=[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]")
        pml_lines.append(f"show sphere, true_center_{i+1}")
        
        # c) real residues selection (no color)
        pml_lines.append(to_pymol_selection(f"true_pocket_{i+1}", target_chain, tp))
    
    pml_lines.append("color green, true_center_*")
    pml_lines.append("")

    # Predicted Pockets processing
    pml_lines.append("# --- TOP N PREDICTED POCKETS ---")
    for i, pp in enumerate(top_n_preds):
        center = resolve_center(pp, target_pdb, target_chain, ca_coords_cache, args.structures)
        
        # d) centers of top n predicted
        pml_lines.append(f"pseudoatom pred_center_{i+1}, pos=[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]")
        pml_lines.append(f"show sphere, pred_center_{i+1}")
        
        # e) residue selections of top n predicted (no coloring)
        pml_lines.append(to_pymol_selection(f"pred_pocket_{i+1}", target_chain, pp["residues"]))

    pml_lines.append("color magenta, pred_center_*")
    pml_lines.append("")
    pml_lines.append("center protein")

    # Write output
    with open(args.output, "w") as f:
        f.write("\n".join(pml_lines) + "\n")
    
    print(f"Generated PyMOL script: {args.output}")
    print(f"To view, run: pymol {args.output}")

if __name__ == "__main__":
    main()