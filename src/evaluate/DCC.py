import numpy as np

def DCC(predicted_centers, true_centers):
    """
    Compute the Distance Center-to-Center (DCC) metric.
    
    Parameters:
    predicted_centers (list of numpy arrays): List of predicted binding site centers (x, y, z).
    true_centers (list of numpy arrays): List of true binding site centers (x, y, z).
    
    Returns:
    min_DCCs (list of floats): List of minimum DCC values for each true center. Length is equal to the number of true centers: len(true_centers) == len(min_DCCs).
    """
    min_DCCs = []
    for true_center in true_centers:
        min_DCC = float('inf')
        for predicted_center in predicted_centers:
            distance = np.linalg.norm(true_center - predicted_center)
            min_DCC = min(min_DCC, distance)
        min_DCCs.append(min_DCC)

    return min_DCCs