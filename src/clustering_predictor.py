import os

import pandas as pd
import numpy as np
import os
import multiprocessing

from sklearn.cluster import DBSCAN
from pocket_pipeline import prediction_pipeline, create_base_parser
from functools import partial

def add_args(parser):
    parser.add_argument("--columns", default="score", help="Comma-separated list of columns to average. Default: score")
    parser.add_argument("--probability_agg", default="mean", choices=["mean", "max"], help="Method to aggregate probabilities across frames. Default: mean")
    parser.add_argument("--weight_type", default="cosine", choices=["cosine", "linear", "none"], help="Type of weighting to apply when averaging scores across frames. Default: cosine")
    parser.add_argument("--eps", type=float, default=2.0, help="DBSCAN eps parameter for clustering pockets. Default: 2.0")
    return parser

def string_list_union(inputs, sep=' '):
    if not inputs:
        return []
    return list(set([val for sublist in inputs for val in sublist.split(sep)]))

def dbscan_cluster(coords, eps=2.0):
    return DBSCAN(eps=eps, min_samples=1, n_jobs=multiprocessing.cpu_count()).fit_predict(coords)

def aggregate_pockets(prediction_df: pd.DataFrame, cluster_func=dbscan_cluster, verbose=0):
    """
    Aggregate pocket predictions across frames for a given protein directory.
    Returns a DataFrame with aggregated pocket information. The cluster function should 
    accept a list of 3D coordinates and return a cluster label for each pocket, which will be saved the 'cluster' column.
    The final output is a DataFrame aggregated by the 'cluster' column.

    Parameters:
    prediction_df (pd.DataFrame): DataFrame containing pocket predictions across frames.
    cluster_func (function): Function that takes in 3D coordinates and returns cluster labels. Default is DBSCAN clustering.

    Returns: pd.DataFrame: DataFrame with aggregated pocket information, including averaged scores and coordinates.
    """

    # Cluster predictions based on spatial proximity
    coords = prediction_df[["center_x", "center_y", "center_z"]].values

    # Using DBSCAN to cluster predictions that are within 2.0 units (assuming Ångströms) of each other
    prediction_df['cluster'] = cluster_func(coords)

    # Aggregate scores by cluster   
    aggregated = prediction_df.groupby('cluster').agg(list).reset_index(drop=True)
    return aggregated

def cosine_weighted_average(scores, n_total_frames, n_cluster_frames):
    """
    Calculate the cosine-weighted average of a list of scores.
    Pockets present only in one frame will have a weight of 1, while pockets present in all frames will have a weight of ~0 (but > 0), with a smooth cosine transition in between.

    Parameters:
    scores (list): List of scores.
    n_total_frames (int): Total number of frames.
    n_cluster_frames (int): Number of frames in the cluster.

    Returns:
    float: Cosine-weighted average
    """
    if len(scores) == 0:
        return 0
    return np.mean(scores) * (1 + np.cos(np.pi * (n_cluster_frames - 1) / n_total_frames)) / 2

def linear_weighted_average(scores, n_total_frames, n_cluster_frames):
    """
    Calculate the linearly-weighted average of a list of scores.
    Pockets present only in one frame will have a weight of 1, while pockets present in all frames will have a weight of 0, with a linear transition in between.

    Parameters:
    scores (list): List of scores.
    n_total_frames (int): Total number of frames.
    n_cluster_frames (int): Number of frames in the cluster.

    Returns:
    float: Linearly-weighted average
    """
    
    if len(scores) == 0:
        return 0
    return np.mean(scores) * (1 - (n_cluster_frames - 1) / n_total_frames)

def no_weight_average(scores, n_total_frames, n_cluster_frames):
    """
    Calculate the simple average of a list of scores without any weighting.

    Parameters:
    scores (list): List of scores.
    n_frames (int): Number of frames (not used in this function).
    
    Returns:
    float: Simple average
    """
    if len(scores) == 0:
        return 0
    return np.mean(scores)

def get_weight_func(weight_type):
    if weight_type == "cosine":
        return cosine_weighted_average
    elif weight_type == "linear":
        return linear_weighted_average
    elif weight_type == "none":
        return no_weight_average
    else:
        raise ValueError(f"Invalid weight type: {weight_type}. Must be one of: cosine, linear, none.")
    
def process_pockets(aggregated_df, weight_func, by_column="score", n_frames=None, verbose=0):
    if by_column not in aggregated_df.columns:
        raise ValueError(f"Column '{by_column}' not found in aggregated DataFrame.")
    
    # Score is the average score across frames for the pocket, weighted by the number of predictions in each cluster
    aggregated_df['n_cluster_frames'] = aggregated_df['frame_file'].apply(lambda x: len(set(x)))
    aggregated_df[by_column] = aggregated_df[[by_column, 'n_cluster_frames']].apply(lambda x: weight_func(x[by_column], n_frames, x['n_cluster_frames']), axis=1)
    for coord in ["center_x", "center_y", "center_z"]:
        aggregated_df[coord] = aggregated_df[coord].apply(lambda x: np.mean(x))
    aggregated_df['residue_ids'] = aggregated_df['residue_ids'].apply(string_list_union)
    aggregated_df = aggregated_df.sort_values(by=by_column, ascending=False)
    aggregated_df['rank'] = range(1, len(aggregated_df) + 1)
    # Just a placeholder really
    if args.probability_agg == "mean":
        aggregated_df['probability'] = aggregated_df['probability'].apply(lambda x: np.mean(x) if isinstance(x, list) else x)
    elif args.probability_agg == "max":
        aggregated_df['probability'] = aggregated_df['probability'].apply(lambda x: np.max(x) if isinstance(x, list) else x)
    aggregated_df['name'] = aggregated_df['rank'].apply(lambda x: f"pocket{x}")
    return aggregated_df.reset_index(drop=True)[['name', 'rank', by_column, 'probability', 'center_x', 'center_y', 'center_z', 'residue_ids']]

if __name__ == "__main__":
    parser = create_base_parser()
    parser = add_args(parser)
    args = parser.parse_args()
    output_path = args.output_dir if args.output_dir else os.path.join(os.path.pardir, "output", f"clustered_{args.weight_type}_weighted_predictions_eps{args.eps}")
    cluster_function = partial(dbscan_cluster, eps=args.eps)
    weight_func = get_weight_func(args.weight_type)
    prediction_pipeline(args, partial(aggregate_pockets, cluster_func=cluster_function), partial(process_pockets, weight_func=weight_func, by_column='score'), output_path=output_path)