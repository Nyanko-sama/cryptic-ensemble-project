
def RRO(cluster_residues, actual_binding_sites):
    """Compute the Residue Recovery Overlap (RRO) metric.
    Args:
        cluster_residues (list of lists): A list where each element is a set of residues predicted by the pLM for a particular cluster.
        actual_binding_sites (list of lists): A list where each element is a set of residues that are known to be part of actual binding sites.
    Returns:
        list: The best computed RRO value for each actual binding site.
    """
    rro = []
    for pocket in actual_binding_sites:
        best_residue_overlap = 0
        # compute RRO for pLM predictions
        for residues in cluster_residues:
            residue_overlap = len(set(pocket).intersection(residues)) / len(pocket)
            if residue_overlap > best_residue_overlap:
                best_residue_overlap = residue_overlap
        rro.append(best_residue_overlap)

    return rro