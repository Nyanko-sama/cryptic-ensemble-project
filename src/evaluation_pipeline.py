import pandas as pd
import numpy as np
import os
import argparse
import json
import sys
import biotite

sys.path.append(os.path.join('..', "utils",))
sys.path.append(os.path.join('..', ''))
sys.path.append(os.path.join('..', 'deeplife_2026', ''))

from utils.data_loading import load_all_predictions, load_evaluation_dataset
from sklearn.metrics import precision_recall_curve, precision_score, recall_score, roc_auc_score, average_precision_score, matthews_corrcoef, f1_score, confusion_matrix
from deeplife_2026.src.DCC import DCC
from deeplife_2026.src.RRO import RRO

def create_parser():
    parser = argparse.ArgumentParser(description="Evaluate performance of a pocket prediction method against the CryptoBench dataset.")
    parser.add_argument("--predictions_json", required=True, help="Path to JSON file containing pocket predictions to evaluate. Must contain columns: name, score, center_x, center_y, center_z, residue_ids, frame_file.")
    parser.add_argument("--crypto_path", default="../data/cryptobench/", help="Path to the cryptobench repository. Default: ../data/cryptobench/")
    parser.add_argument("--eval_dataset_path", default="../data/test_eval_dataset_auth_labels.csv", help="Path to the evaluation dataset CSV file. Default: ../data/deeplife_2026/test.csv")
    parser.add_argument("--DCC_threshold", type=float, default=12, help="Threshold for DCC metric. Default: 5.0")
    parser.add_argument("--k", type=int, default=0, help="Number of top predictions to consider. Default: 0")
    return parser


def create_graphs(grouped_predictions_df, eval_dataset):
    pass

def compute_single_metrics(predictions, labels):
    precision = precision_score(labels, predictions)
    recall = recall_score(labels, predictions)
    f1 = f1_score(labels, predictions)
    mcc = matthews_corrcoef(labels, predictions)
    auc_roc = roc_auc_score(labels, predictions)
    auc_pr = average_precision_score(labels, predictions)

    return {
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'mcc': mcc,
        'auc_roc': auc_roc,
        'auc_pr': auc_pr,
    }

def compute_metrics(grouped_predictions_df : pd.DataFrame, eval_dataset : pd.DataFrame, DCC_threshold, k=0):
    # Merge predictions with evaluation dataset on pdb_id and chain_id
    # NOTE: Chain ids will be added later, for now use only pdb_id for merging and ensure consistency in the dataset
    grouped_predictions_df['centers'] = grouped_predictions_df.apply(lambda row: np.array([row['center_x'], row['center_y'], row['center_z']]).T, axis=1)
    merged_df = grouped_predictions_df.merge(eval_dataset, left_on=['pdb_id'], right_on=['pdb_id'], how='inner', suffixes=('_pred', '_eval'))
    if grouped_predictions_df.shape[0] != merged_df.shape[0]:
        print("Warning: Some predictions could not be matched with the evaluation dataset. Check pdb_id and chain_id consistency.")

    # Compute DCC and RRO for each prediction
    merged_df['DCC'] = merged_df.apply(lambda row: DCC(row['centers_pred'][:len(row['centers_pred']) + k], row['centers_eval']), axis=1)
    merged_df['RRO'] = merged_df.apply(lambda row: RRO(row['residue_ids'][:len(row['residue_ids']) + k], row['binding_residue_ids']), axis=1)

    dccs = np.array([val for lst in merged_df['DCC'].to_list() for val in lst])
    rros = np.array([val for lst in merged_df['RRO'].to_list() for val in lst])

    prefix = f'TopN' if k == 0 else f'TopN+{k}'
    metrics = {
        f'{prefix}_mean_DCC': dccs.mean(),
        f'{prefix}_acc' : np.sum(dccs < DCC_threshold)/dccs.shape[0],
        f'{prefix}_median_DCC': np.median(dccs),
        f'{prefix}_mean_RRO': rros.mean(),
        f'{prefix}_median_RRO': np.median(rros),
    }
    return metrics
    
def save_metrics(metrics, output_path):
    with open(output_path, 'w') as f:
        json.dump(metrics, f, indent=4)

def nice_print_metrics(metrics):
    for metric, value in metrics.items():
        print(f"{metric}: {value}")

def evaluation_pipeline(predictions_json, eval_dataset_path, DCC_threshold=4.0, k=0):
    # Load predictions and target dataset
    grouped_predictions_df = load_all_predictions(predictions_json)
    print(grouped_predictions_df.head(1))
    # assumes the dataset already passed through prepare_dataset.py and contains 'binding_residue_ids' and 'residue_ids' columns
    grouped_eval_dataset = load_evaluation_dataset(eval_dataset_path)
    print(grouped_eval_dataset.head(1))

    # Compute metrics
    metrics = compute_metrics(grouped_predictions_df, grouped_eval_dataset, DCC_threshold=DCC_threshold, k=k)
    nice_print_metrics(metrics)
    return metrics



if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    metrics = evaluation_pipeline(args.predictions_json, args.eval_dataset_path, DCC_threshold=args.DCC_threshold, k=args.k)
    nice_print_metrics(metrics)
    args.output_path = os.path.join(os.path.dirname(args.predictions_json), f"evaluation_metrics_topn_plus{args.k}_{args.DCC_threshold}.json")
    save_metrics(metrics, args.output_path)