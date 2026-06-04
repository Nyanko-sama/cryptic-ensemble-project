def FPI(pred_residues_list, true_residues_list):
    """
    For each true pocket, count how many predicted pockets (MAX mode) have any
    residue overlap with it. Returns those counts for true pockets that have at
    least one overlapping predicted pocket (others are skipped).

    Aggregate over all proteins: mean of returned counts.
    Example: FPI=2 means each covered true pocket is overlapped by 2 predictions
    on average.

    Parameters
    ----------
    pred_residues_list  : list of collections  (all predicted pockets, MAX mode)
    true_residues_list  : list of collections  (true pockets for this protein)

    Returns
    -------
    list of int  — one count per true pocket with ≥1 overlap (may be empty)
    """
    counts = []
    for true_res in true_residues_list:
        true_set = set(true_res)
        count = sum(1 for pred_res in pred_residues_list if true_set.intersection(pred_res))
        if count > 0:
            counts.append(count)
    return counts
