import pandas as pd
import numpy as np
import os
import argparse
import sys
import multiprocessing

from sklearn.cluster import DBSCAN
from pocket_pipeline import prediction_pipeline, create_base_parser
from functools import partial

def add_args(parser):
    parser.add_argument("--columns", default="score", help="Comma-separated list of columns to average. Default: score")

    return parser
def string_list_union(list, sep=';'):
    if not list:
        return []
    return list(set([val for sublist in list for val in sublist.split(sep)]))

def dbscan_cluster(coords, eps=2.0):
    return DBSCAN(eps=eps, min_samples=1, n_jobs=multiprocessing.cpu_count()).fit_predict(coords)

def aggregate_pockets(prediction_df: pd.DataFrame, cluster_function=dbscan_cluster):
    """
    Aggregate pocket predictions across frames for a given protein directory.
    Returns a DataFrame with aggregated pocket information. The cluster function should 
    accept a list of 3D coordinates and return a cluster label for each pocket, which will be saved the 'cluster' column.
    The final output is a DataFrame aggregated by the 'cluster' column.

    Parameters:
    prediction_df (pd.DataFrame): DataFrame containing pocket predictions across frames.
    cluster_function (function): Function that takes in 3D coordinates and returns cluster labels. Default is DBSCAN clustering.

    Returns: pd.DataFrame: DataFrame with aggregated pocket information, including averaged scores and coordinates.
    """
    

    # Cluster predictions based on spatial proximity
    coords = prediction_df[["center_x", "center_y", "center_z"]].values

    # Using DBSCAN to cluster predictions that are within 2.0 units (assuming Ångströms) of each other
    prediction_df['cluster'] = cluster_function(coords, eps=2.0)

    # Aggregate scores by cluster   
    aggregated = prediction_df.groupby('cluster').agg(list).reset_index(drop=True)
    return aggregated

def cosine_weighted_average(scores, n_frames):
    """
    Calculate the cosine-weighted average of a list of scores.
    Pockets present only in one frame will have a weight of 1, while pockets present in all frames will have a weight of ~0 (but > 0), with a smooth cosine transition in between.

    Parameters:
    scores (list): List of scores.
    n_frames (int): Number of frames.
    
    Returns:
    float: Cosine-weighted average
    """
    if not scores:
        return 0
    return np.mean(scores) * (1 + np.cos(np.pi * (len(scores) - 1) / n_frames)) / 2

def process_pockets(aggregated_df, by_column="score", n_frames=None):
    if by_column not in aggregated_df.columns:
        raise ValueError(f"Column '{by_column}' not found in aggregated DataFrame.")
    
    # Score is the average score across frames for the pocket, weighted by the number of predictions in each cluster
    aggregated_df[by_column] = aggregated_df[by_column].apply(lambda x: cosine_weighted_average(x, n_frames))
    for coord in ["center_x", "center_y", "center_z"]:
        aggregated_df[coord] = aggregated_df[coord].apply(lambda x: np.mean(x))
    aggregated_df['residue_ids'] = aggregated_df['residue_ids'].apply(string_list_union)
    aggregated_df = aggregated_df.drop(columns=['cluster'])
    aggregated_df = aggregated_df.sort_values(by=by_column, ascending=False)
    aggregated_df['rank'] = range(1, len(aggregated_df) + 1)
    
    return aggregated_df

if __name__ == "__main__":
    parser = create_base_parser()
    parser = add_args(parser)
    args = parser.parse_args()
    output_path = args.output_path if args.output_path else args.base_dir
    prediction_pipeline(args, partial(aggregate_pockets, cluster_function=dbscan_cluster), partial(process_pockets, by_column='score'), output_path=output_path)