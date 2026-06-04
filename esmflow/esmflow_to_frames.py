import os
import argparse

def esmfold_to_pdb_text(esmfold_pdb_path: str, output_dir: str) -> list[str]:
    """Convert a multi-model ESMFold PDB into separate single-frame PDB files using text parsing."""
    os.makedirs(output_dir, exist_ok=True)
    proteins_paths = []
    
    with open(esmfold_pdb_path, 'r') as f:
        lines = f.readlines()

    current_frame_lines = []
    frame_idx = 0
    in_model = False
    
    for line in lines:
        if line.startswith("MODEL"):
            in_model = True
            current_frame_lines = []  # Clear previous lines for the new model
            continue  # Skip writing the MODEL line itself
            
        if line.startswith("ENDMDL"):
            # Save accumulated lines to a new PDB file
            out_path = os.path.join(output_dir, f"frame_{frame_idx:05d}.pdb")
            with open(out_path, 'w') as out_f:
                out_f.writelines(current_frame_lines)
            
            proteins_paths.append(out_path)
            frame_idx += 1
            in_model = False
            continue
            
        if in_model:
            current_frame_lines.append(line)

    # Fallback in case the very last model is missing an ENDMDL tag
    if in_model and current_frame_lines:
        out_path = os.path.join(output_dir, f"frame_{frame_idx:05d}.pdb")
        with open(out_path, 'w') as out_f:
            out_f.writelines(current_frame_lines)
        proteins_paths.append(out_path)

    return proteins_paths


def main():
    parser = argparse.ArgumentParser(description="Convert ESMFold multi-model PDB to single-frame PDBs.")
    parser.add_argument("--esmflow_dir", type=str, required=True, help="Path to the ESMFold multi-model PDB file.")
    
    args = parser.parse_args()

    for folder in os.listdir(args.esmflow_dir):
        path = os.path.join(args.esmflow_dir, folder)
        if os.path.isdir(path):
            esmfold_pdb_path = os.path.join(path, f"{folder}.pdb")
            
            # Ensure the PDB file actually exists before trying to parse it
            if not os.path.exists(esmfold_pdb_path):
                print(f"Warning: {esmfold_pdb_path} not found. Skipping.")
                continue

            frame_paths = esmfold_to_pdb_text(esmfold_pdb_path, path)
            print(f"Generated {len(frame_paths)} frame PDB files in {path}")


if __name__ == "__main__":
    main()