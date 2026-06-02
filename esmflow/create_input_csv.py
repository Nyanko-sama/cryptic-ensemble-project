import json
import csv
import os

def create_csv_from_json(json_path, csv_path):
    """
    Reads a JSON file with sequences and creates a CSV file.

    The JSON file is expected to have a dictionary where keys are in the
    format '<pdb_id>_<chain>_<model>.pdb' and values are sequences.

    The output CSV will have two columns: 'name' (the pdb_id) and 'seqres' (the sequence).
    """
    with open(json_path, 'r') as f:
        data = json.load(f)

    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'seqres'])
        for key, value in data.items():
            pdb_id = key.split('_')[0]
            writer.writerow([pdb_id, value['sequence']])

if __name__ == '__main__':
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct absolute paths
    json_file = os.path.join(script_dir, '..', 'data', 'cryptobench', 'cryptobench-dataset', 'sequences.json')
    csv_file = os.path.join(script_dir, 'input.csv')
    
    create_csv_from_json(json_file, csv_file)
    print(f"Successfully created {csv_file} from {json_file}")
