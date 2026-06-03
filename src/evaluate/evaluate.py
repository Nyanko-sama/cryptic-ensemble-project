#!/usr/bin/env python3
"""
evaluate.py — DeepLife 2026 Submission Evaluator

Usage:
    python src/evaluate.py \\
        --predictions src/predictions.json \\
        --ground-truth data/test.csv \\
        --structures data/structures/ \\
        [--dcc-threshold 12.0] \\
        [--topn-k 2] \\
        [--output results.json]

Evaluation modes
----------------
Top-N    : candidates = top-N ranked pockets, N = number of true pockets for that protein
Top-N+K  : candidates = top-(N+K) ranked pockets  (K set by --topn-k, default 2)
MAX      : candidates = all submitted pockets (oracle upper bound)

Metrics
-------
DCC success rate : fraction of true pockets with assigned DCC < threshold.
                   1:1 assignment (Hungarian); denominator = all GT pockets.
FPI              : fraction of proteins where ≥1 true pocket has DCC < threshold.
                   Denominator = all GT proteins (unsubmitted = miss).
RRO              : mean residue-recovery overlap over pockets with overlap > 0.
                   1:1 assignment (Hungarian).

If 'center' is absent from a pocket, the Cα centroid is computed from the
pocket's residue list using the CIF structures (--structures required).
"""

import argparse
import json
import sys
from pathlib import Path
import numpy as np

_SRC = Path(__file__).parent
sys.path.insert(0, str(_SRC))
from DCC import DCC as compute_dcc
from RRO import RRO as compute_rro
from FPI import FPI as compute_fpi


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_ground_truth(csv_path):
    """
    Returns dict: (pdb_id_lower, chain) -> list[list[str]]
    Each inner list is the auth_seq_id residues for one pocket as strings.
    Format: PDB_ID;Chain;Ligands;Residue_Indices  (semicolon-delimited)
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
    Converts ["B:220", "A:-3", "A:123B"] -> ["220", "-3", "123B"].
    Strips chain prefix; preserves insertion codes.
    """
    return [r.split(":", 1)[1] for r in residue_list]


# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def _load_ca_coords(cif_path, chain_id):
    """
    Returns dict {residue_key: np.array([x,y,z])} for all Cα in chain.
    Keys are strings like "220" or "60A" (res_id + ins_code).
    """
    import biotite.structure.io.pdbx as pdbx
    from biotite.structure.io.pdbx import get_structure
    from biotite.structure import get_residues

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

    assert len(residue_ids) == len(coords)
    return coords


def _centroid(ca_coords, residue_keys):
    pts = [ca_coords[r] for r in residue_keys if r in ca_coords]
    if not pts:
        raise ValueError(f"No Cα atoms found for residues {residue_keys}")
    return np.mean(pts, axis=0)


def _get_ca_cache(pdb_id, chain, ca_coords_cache, structures_dir):
    key = (pdb_id.lower(), chain)
    if key not in ca_coords_cache:
        if structures_dir is None:
            raise ValueError(
                f"--structures is required to compute centers for {pdb_id}/{chain}"
            )
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
    return _centroid(ca, parse_residues(pocket["residues"]))


def gt_center(residue_keys, pdb_id, chain, ca_coords_cache, structures_dir):
    ca = _get_ca_cache(pdb_id, chain, ca_coords_cache, structures_dir)
    return _centroid(ca, residue_keys)


# ---------------------------------------------------------------------------
# Per-protein evaluation
# ---------------------------------------------------------------------------

def evaluate_protein(ranked_pockets, true_pockets, pdb_id, chain,
                     structures_dir, ca_coords_cache, topn_k):
    """
    Returns:
        dccs_n, dccs_nk, dccs_max  — DCC per true pocket (closest-neighbor)
        rro_max                    — RRO per true pocket over all candidates
        fpi_counts                 — overlap counts per covered true pocket (MAX)
    """
    n_true = len(true_pockets)
    n_pred = len(ranked_pockets)

    pred_centers = []
    pred_residues = []
    for p in ranked_pockets:
        pred_centers.append(resolve_center(p, pdb_id, chain, ca_coords_cache, structures_dir))
        pred_residues.append(parse_residues(p["residues"]))

    true_centers = [
        gt_center(tp, pdb_id, chain, ca_coords_cache, structures_dir)
        for tp in true_pockets
    ]
    true_residues = [list(tp) for tp in true_pockets]

    n_for_n   = min(n_true, n_pred)
    n_for_nk  = min(n_true + topn_k, n_pred)
    n_for_max = n_pred

    dccs_n   = compute_dcc(pred_centers[:n_for_n],   true_centers)
    dccs_nk  = compute_dcc(pred_centers[:n_for_nk],  true_centers)
    dccs_max = compute_dcc(pred_centers[:n_for_max], true_centers)

    rro_max    = compute_rro(pred_residues[:n_for_max], true_residues)
    fpi_counts = compute_fpi(pred_residues[:n_for_max], true_residues)

    return dccs_n, dccs_nk, dccs_max, rro_max, fpi_counts


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Evaluate DeepLife 2026 cryptic pocket predictions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--predictions", required=True,
                   help="Path to predictions.json")
    p.add_argument("--ground-truth", required=True,
                   help="Path to test.csv (semicolon-delimited)")
    p.add_argument("--structures", default=None,
                   help="Directory with .cif files; required when 'center' is "
                        "absent from any pocket in the submission")
    p.add_argument("--dcc-threshold", type=float, default=12.0,
                   help="Primary DCC success threshold in Å")
    p.add_argument("--topn-k", type=int, default=2,
                   help="K offset for Top-N+K mode")
    p.add_argument("--output", default=None,
                   help="Write full results to this JSON file")
    return p.parse_args()


def main():
    args = parse_args()
    thr = args.dcc_threshold
    K   = args.topn_k

    gt = load_ground_truth(args.ground_truth)
    submission = load_submission(args.predictions)

    total_pockets  = sum(len(v) for v in gt.values())
    total_proteins = len(gt)

    ca_coords_cache = {}
    all_dccs_n,  all_dccs_nk,  all_dccs_max  = [], [], []
    all_rro_max = []
    all_fpi_counts = []
    skipped = []

    for pred in submission["predictions"]:
        pdb_id = pred["pdb_id"].lower()
        chain  = pred["chain"]
        key    = (pdb_id, chain)

        if key not in gt:
            skipped.append(f"{pdb_id}/{chain} — not in ground truth")
            continue

        ranked = sorted(pred["ranked_pockets"], key=lambda p: p["rank"])

        try:
            dccs_n, dccs_nk, dccs_max, rro_max, fpi_counts = evaluate_protein(
                ranked, gt[key], pdb_id, chain,
                args.structures, ca_coords_cache, K,
            )
        except Exception as e:
            skipped.append(f"{pdb_id}/{chain} — {e}")
            print(f"  WARNING: {pdb_id}/{chain}: {e}", file=sys.stderr)
            continue

        all_dccs_n.extend(dccs_n)
        all_dccs_nk.extend(dccs_nk)
        all_dccs_max.extend(dccs_max)
        all_rro_max.extend(v for v in rro_max if v > 0)
        all_fpi_counts.extend(fpi_counts)

    def success_rate(dccs):
        return float(np.sum(np.array(dccs) < thr) / total_pockets) if total_pockets else 0.0

    def _mean(vals):
        return float(np.mean(vals))   if vals else float("nan")

    def _median(vals):
        return float(np.median(vals)) if vals else float("nan")

    dcc_n   = success_rate(all_dccs_n)
    dcc_nk  = success_rate(all_dccs_nk)
    dcc_max = success_rate(all_dccs_max)
    rro_mean   = _mean(all_rro_max)
    rro_median = _median(all_rro_max)
    fpi_mean   = _mean(all_fpi_counts)

    n_evaluated = len(submission["predictions"]) - len(skipped)

    print(f"\nnumber_of_pockets (GT total): {total_pockets}")
    print(f"number_of_proteins (GT total): {total_proteins}")
    print(f"proteins evaluated: {n_evaluated}  |  skipped: {len(skipped)}")
    print()
    col = f"{'Mode':<8}  {'DCC':>8}"
    sep = "-" * len(col)
    print(col)
    print(sep)
    print(f"{'MAX':<8}  {dcc_max:>8.4f}")
    print(f"{'N+'+str(K):<8}  {dcc_nk:>8.4f}")
    print(f"{'N':<8}  {dcc_n:>8.4f}")
    print(f"\nRRO (MAX):  mean={rro_mean:.4f},  median={rro_median:.4f}")
    print(f"FPI (MAX):  mean={fpi_mean:.4f}  (n={len(all_fpi_counts)} covered pockets)")

    if skipped:
        print(f"\nSkipped ({len(skipped)}):")
        for s in skipped:
            print(f"  - {s}")

    output_data = {
        "metadata": submission.get("metadata", {}),
        "dcc_threshold": thr,
        "topn_k": K,
        "total_gt_pockets": total_pockets,
        "total_gt_proteins": total_proteins,
        "n_proteins_evaluated": n_evaluated,
        "aggregate": {
            "dcc_max":         {"success_rate": dcc_max},
            f"dcc_n_plus_{K}": {"success_rate": dcc_nk},
            "dcc_n":           {"success_rate": dcc_n},
            "rro_max":         {"mean": rro_mean, "median": rro_median},
            "fpi_max":         {"mean": fpi_mean, "n_covered_pockets": len(all_fpi_counts)},
        },
    }

    if args.output:
        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults written to {args.output}")

    return output_data


if __name__ == "__main__":
    main()
