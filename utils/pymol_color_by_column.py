from __future__ import annotations

import argparse
import csv
import os
import sys

from pymol import cmd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Color a PyMOL object by avg_probability values from residue_averages.csv."
    )
    parser.add_argument(
        "csv_file",
        help="Path to residue_averages.csv containing avg_probability.",
    )
    parser.add_argument(
        "object_name",
        nargs="?",
        default=None,
        help="PyMOL object name to color. If omitted, the first loaded object is used.",
    )
    parser.add_argument(
        "--column",
        default="avg_probability",
        help="Column name in CSV to use for coloring (default: avg_probability).",
    )
    parser.add_argument(
        "--selection",
        default="all",
        help="PyMOL selection to color (default: all).",
    )
    parser.add_argument(
        "--palette",
        default="blue_white_red",
        help="PyMOL color palette or ramp name to use (default: blue_white_red).",
    )
    return parser.parse_args()


def find_default_object() -> str:
    objects = cmd.get_names("objects")
    if not objects:
        raise RuntimeError("No PyMOL objects are loaded.")
    return objects[0]


def read_residue_averages(csv_path: str, column: str) -> list[tuple[str, str, float]]:
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    residues: list[tuple[str, str, float]] = []
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        headers = [h.strip() for h in reader.fieldnames or []]
        if column not in headers:
            raise ValueError(
                f"Column '{column}' not found in CSV. Available columns: {', '.join(headers)}"
            )
        for row in reader:
            chain = row.get("chain", "").strip()
            residue_label = row.get("residue_label", "").strip()
            if not chain or not residue_label:
                continue
            value_str = row.get(column, "").strip()
            if not value_str:
                continue
            try:
                value = float(value_str)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid numeric value '{value_str}' for {column} on chain {chain}, residue {residue_label}"
                ) from exc
            residues.append((chain, residue_label, value))
    if not residues:
        raise ValueError(f"No valid residue entries found in {csv_path}.")
    return residues


def set_residue_b_factors(object_name: str, residues: list[tuple[str, str, float]]) -> tuple[float, float]:
    values = []
    for chain, residue_label, value in residues:
        selection = f"{object_name} and chain {chain} and resi {residue_label}"
        if cmd.count_atoms(selection) == 0:
            continue
        cmd.alter(selection, f"b={value}")
        values.append(value)

    if not values:
        raise RuntimeError("No matching residues found in the object for the provided CSV entries.")

    return min(values), max(values)


def color_by_b(object_name: str, selection: str, palette: str, min_val: float, max_val: float) -> None:
    full_selection = f"{object_name} and ({selection})"
    cmd.spectrum("b", palette, full_selection, min_val, max_val)
    cmd.cartoon("putty", full_selection)
    cmd.set("cartoon_transparency", 0.0, full_selection)
    cmd.set("cartoon_putty_scale_min", 0.5, full_selection)
    cmd.set("cartoon_putty_scale_max", 2.0, full_selection)


def pymol_color_by_column(csv_file: str, object_name: str = None, column: str = "avg_probability", selection: str = "all", palette: str = "blue_white_red") -> None:
    """PyMOL command wrapper for coloring by a CSV column."""
    if object_name is None:
        object_name = find_default_object()

    residues = read_residue_averages(csv_file, column)
    min_val, max_val = set_residue_b_factors(object_name, residues)
    color_by_b(object_name, selection, palette, min_val, max_val)
    print(
        f"Colored object '{object_name}' by {column} using {len(residues)} residues."
        f" Range: {min_val:.6g} to {max_val:.6g}"
    )


cmd.extend("pymol_color_by_column", pymol_color_by_column)


def main() -> None:
    args = parse_args()
    if args.object_name is None:
        args.object_name = find_default_object()

    residues = read_residue_averages(args.csv_file, args.column)
    min_val, max_val = set_residue_b_factors(args.object_name, residues)
    color_by_b(args.object_name, args.selection, args.palette, min_val, max_val)
    print(
        f"Colored object '{args.object_name}' by {args.column} using {len(residues)} residues."
        f" Range: {min_val:.6g} to {max_val:.6g}"
    )


if __name__ == "__main__":
    main()
