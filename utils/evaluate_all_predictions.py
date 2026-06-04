#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import os


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
    parser.add_argument(
        '--jobs',
        type=int,
        default=max(1, (os.cpu_count() or 2) // 2),
        help='Number of parallel evaluation jobs to run. Defaults to half the CPUs.',
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Re-run evaluations even if output files already exist.',
    )
    return parser.parse_args()


def main():
    args = parse_args()
    predictions_files = find_all_predictions(args.output_root)

    if not predictions_files:
        print(f"No all_predictions.json files found under {args.output_root}")
        raise SystemExit(1)

    # Collect tasks and run in a thread pool so subprocesses run concurrently.
    futures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for predictions_path in predictions_files:
            parent_dir = predictions_path.parent
            train_output = parent_dir / 'all_predictions.eval_train.json'
            test_output = parent_dir / 'all_predictions.eval_test.json'

            # If both evaluation outputs exist and not forcing, skip this predictions file.
            if not args.force and train_output.exists() and test_output.exists():
                print(f"Skipping evaluation for {predictions_path} — both train and test eval files exist")
                continue

            # Schedule train evaluation if forcing or missing
            if args.force or not train_output.exists():
                futures.append(pool.submit(
                    run_evaluation,
                    args.evaluate_script,
                    predictions_path,
                    args.train_ground_truth,
                    args.structures,
                    train_output,
                ))
            else:
                print(f"Train evaluation already exists for {predictions_path}, skipping train")

            # Schedule test evaluation if forcing or missing
            if args.force or not test_output.exists():
                futures.append(pool.submit(
                    run_evaluation,
                    args.evaluate_script,
                    predictions_path,
                    args.test_ground_truth,
                    args.structures,
                    test_output,
                ))
            else:
                print(f"Test evaluation already exists for {predictions_path}, skipping test")

        # Wait for all scheduled evaluations to finish and propagate exceptions
        for fut in as_completed(futures):
            try:
                fut.result()
            except SystemExit as e:
                print(f"Evaluation subprocess exited with code {e.code}")
                raise
            except Exception:
                print("An evaluation job failed:")
                raise

    print('Evaluation complete.')


if __name__ == '__main__':
    main()
