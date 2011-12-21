from scipy import *
from numpy import *

def compute_multi_iso_credible_intervals(multi_iso_samples,
                                         confidence_level=0.95):
    """
    Compute multiple isoforms credible intervals for a set of NxM matrix.
    """
    credible_intervals = []
    num_samples, num_isoforms = shape(multi_iso_samples)

    for iso_num in range(num_isoforms):
        ci = compute_credible_intervals(multi_iso_samples[:, iso_num],
                                        confidence_level=confidence_level)
        credible_intervals.append(ci)
        
    return credible_intervals
        
