import argparse
import csv
import json
from pathlib import Path


def load_chain_lookup(train_csv: Path, test_csv: Path) -> dict:
    chain_lookup = {}
    for csv_file in (train_csv, test_csv):
        if not csv_file.exists():
            continue

        with csv_file.open(newline='') as f:
            reader = csv.reader(f, delimiter=';')
            for row in reader:
                if len(row) < 2:
                    continue
                pdb_id = row[0].strip()
                chain = row[1].strip()
                if pdb_id:
                    if pdb_id in chain_lookup and chain_lookup[pdb_id] != chain:
                        print(f"Warning: Conflicting chain information for {pdb_id} in {csv_file}. Existing: {chain_lookup[pdb_id]}, New: {chain}. Using existing value.")
                    
                    chain_lookup[pdb_id] = chain

    return chain_lookup


def main(data_dir: Path, output_path: Path) -> None:
    train_csv = data_dir / 'train.csv'
    test_csv = data_dir / 'test.csv'
    chain_lookup = load_chain_lookup(train_csv, test_csv)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        json.dump(chain_lookup, f, indent=2)

    print(f'Wrote chain lookup for {len(chain_lookup)} proteins to {output_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract chain information from train/test CSV files.')
    parser.add_argument(
        '--data-dir',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'data',
        help='Directory containing train.csv and test.csv (default: ../data)'
    )
    parser.add_argument(
        '--output-file',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'data' / 'chain_lookup.json',
        help='Output JSON file path (default: data/chain_lookup.json)'
    )
    args = parser.parse_args()
    main(args.data_dir, args.output_file)
