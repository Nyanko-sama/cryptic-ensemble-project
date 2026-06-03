#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def find_all_predictions(output_root: Path):
    return sorted(output_root.rglob('all_predictions.json'))


def run_evaluation(evaluate_script: Path, predictions_path: Path, ground_truth_path: Path, structures_dir: Path, output_path: Path):
    cmd = [
        'python',
        str(evaluate_script),
        '--predictions', str(predictions_path),
        '--ground-truth', str(ground_truth_path),
        '--structures', str(structures_dir),
        '--output', str(output_path),
    ]
    print(f"Running evaluation for {predictions_path} against {ground_truth_path}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: Evaluation failed for {predictions_path} with ground truth {ground_truth_path}")
        print(result.stdout)
        print(result.stderr)
        raise SystemExit(result.returncode)
    print(result.stdout)
    if result.stderr:
        print(result.stderr)
    print(f"Wrote evaluation results to {output_path}\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description='Evaluate all DeepLife prediction outputs under an output directory for train and test datasets.'
    )
    parser.add_argument(
        '--output-root',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'output',
        help='Root output folder to scan for all_predictions.json files.',
    )
    parser.add_argument(
        '--structures',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'data' / 'cryptobench' / 'cryptobench-dataset' / 'auxiliary-data' / 'cif-files',
        help='Directory containing reference CIF structure files.',
    )
    parser.add_argument(
        '--evaluate-script',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'deeplife_2026' / 'src' / 'evaluate.py',
        help='Path to the DeepLife evaluation script.',
    )
    parser.add_argument(
        '--train-ground-truth',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'data' / 'train.csv',
        help='Path to train ground truth CSV.',
    )
    parser.add_argument(
        '--test-ground-truth',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'data' / 'test.csv',
        help='Path to test ground truth CSV.',
    )
    return parser.parse_args()


def main():
    args = parse_args()
    predictions_files = find_all_predictions(args.output_root)

    if not predictions_files:
        print(f"No all_predictions.json files found under {args.output_root}")
        raise SystemExit(1)

    for predictions_path in predictions_files:
        parent_dir = predictions_path.parent
        train_output = parent_dir / 'all_predictions.eval_train.json'
        test_output = parent_dir / 'all_predictions.eval_test.json'

        run_evaluation(
            args.evaluate_script,
            predictions_path,
            args.train_ground_truth,
            args.structures,
            train_output,
        )
        run_evaluation(
            args.evaluate_script,
            predictions_path,
            args.test_ground_truth,
            args.structures,
            test_output,
        )

    print('Evaluation complete.')


if __name__ == '__main__':
    main()
