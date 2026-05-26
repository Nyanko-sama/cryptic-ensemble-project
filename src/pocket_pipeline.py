
import argparse
import os
import json
from turtle import pd
from utils.data_loading import gather_predictions, get_protein_dirs

def create_base_parser():   
    parser = argparse.ArgumentParser(description="Run the full prediction pipeline for cryptic pocket detection.")
    parser.add_argument("--protein", required=True, help="Protein name or directory containing protein data.")
    parser.add_argument("--base_dir", help="Base directory where protein folders are located. Defaults to ../data/p2rank_preds relative to script.")
    parser.add_argument("--output_path", default="final_predictions.csv", help="Path to save the final predictions CSV file. Default: final_predictions.csv")
    parser.add_argument("--recursive", action="store_true", help="Recursively search for protein directories under the base directory.")
    return parser

def prediction_pipeline(args, aggregation_func, process_func, output_path):
    all_preds = {}
    for protein_dir in get_protein_dirs(args.protein, args.base_dir, recursive=args.recursive):
        predictions_df = gather_predictions(protein_dir)
        if predictions_df.empty:
            raise ValueError(f"No predictions found in protein directory: {protein_dir}")
        n_frames = predictions_df.shape[0]
        aggregation_out = aggregation_func(protein_dir)
        final_pred_df = process_func(aggregation_out, n_frames=n_frames)
        all_preds[os.path.basename(protein_dir)] = final_pred_df.to_dict(orient='records')
        save_predictions(final_pred_df, os.path.join(output_path, os.path.basename(protein_dir) + "_final_predictions.csv"))
    with open(os.path.join(output_path, "all_predictions.json"), 'w') as f:
        json.dump(all_preds, f, indent=4)

def save_predictions(prediction_df, output_path):
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    prediction_df.to_csv(output_path, index=False)