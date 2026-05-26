
import argparse
import os
import sys
import json

sys.path.append(os.path.join('..', 'utils'))

from data_loading import gather_predictions, get_protein_dirs

def create_base_parser():   
    parser = argparse.ArgumentParser(description="Run the full prediction pipeline for cryptic pocket detection.")
    parser.add_argument("--protein", default=None, help="Protein name or directory containing protein data. If None, will process all proteins in the base directory.")
    parser.add_argument("--base_dir", default="../data/p2rank_preds", help="Base directory where protein folders are located.")
    parser.add_argument("--output_path", default=None, help="Path to save the final predictions CSV file.")
    return parser

def prediction_pipeline(args, aggregation_func, process_func, output_path):
    all_preds = {}
    for protein_dir in get_protein_dirs(args.protein, args.base_dir, recursive=args.protein is None):
        predictions_df, n_frames = gather_predictions(protein_dir)
        if predictions_df.empty:
            print(f"No valid predictions found in protein directory: {protein_dir}, skipping.")
            continue
        aggregation_out = aggregation_func(predictions_df)
        final_pred_df = process_func(aggregation_out, n_frames=n_frames)
        all_preds[os.path.basename(protein_dir)] = final_pred_df.to_dict(orient='records')
        save_predictions(final_pred_df, os.path.join(output_path, os.path.basename(protein_dir) + "_final_predictions.csv"))
    save_all_predictions(all_preds, os.path.join(output_path, "all_predictions.json"))

def save_all_predictions(all_preds, output_path):
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    with open(output_path, 'w') as f:
        json.dump(all_preds, f, indent=4)

def save_predictions(prediction_df, output_path):
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    prediction_df.to_csv(output_path, index=False)