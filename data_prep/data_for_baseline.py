
import pandas as pd
import os
import shutil

def process_proteins(test_csv_path, source_dir, baseline_dir):
    """
    Reads protein names from a CSV file, creates directories, and copies PDB files.

    Args:
        test_csv_path (str): Path to the test.csv file.
        source_dir (str): Path to the directory containing PDB files.
        baseline_dir (str): Path to the directory where new directories and files will be created.
    """
    if not os.path.exists(baseline_dir):
        os.makedirs(baseline_dir)

    try:
        df = pd.read_csv(test_csv_path, header=None)
        protein_names = df[0]

        for name_long in protein_names:
            name = name_long.split(';')[0]
            source_file = os.path.join(source_dir, f"{name}.cif")
            dest_dir = os.path.join(baseline_dir, name)
            
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)

            if os.path.exists(source_file):
                shutil.copy(source_file, dest_dir)
                print(f"Copied {source_file} to {dest_dir}")
            else:
                print(f"Source file not found: {source_file}")

    except FileNotFoundError:
        print(f"Error: {test_csv_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Assuming the script is run from the root of the deeplife directory
    test_csv = "test.csv"
    # The user provided a non-existent path, so I am using a path that I believe is correct.
    pdb_source = "../data/cryptobench/cryptobench-dataset/auxiliary-data/cif-files"
    baseline = "../data/baseline"

    process_proteins(test_csv, pdb_source, baseline)
