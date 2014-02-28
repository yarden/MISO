##
## Proposal functions for Metropolis-Hastings
##
## Yarden Katz <yarden@mit.edu>
##
cimport numpy as np
import numpy as np

cimport stat_helpers
cimport matrix_utils
cimport math_utils

def propose_norm_drift_psi_alpha(np.ndarray[double, ndim=1] alpha_vector,
                                 np.ndarray[double, ndim=2] sigma_mat):
    """
    Propose a new alpha vector (parameters to Beta distribution)
    and Psi vector (drawn from a Beta with those alpha vector)
    given a sigma matrix. The new alpha vector is drawn
    from a logistic normal distribution, i.e. by drawing
    a vector from a multivariate normal distribution and then
    taking an inverse logit transform (logit^-1()) to get
    values on the simplex.

    Returns: psi vector (1d), alpha vector (1d)

    if k is number of isoforms, psi vector is of length k
    and alpha vector is of length k-1
    """
    cdef np.ndarray[double, ndim=1] new_psi_vector = \
      np.empty(alpha_vector.shape[0] + 1, dtype=float)
    cdef np.ndarray[double, ndim=1] new_alpha = \
      np.empty(alpha_vector.shape[0], dtype=float)
    if alpha_vector.shape[0] == 1:
        pass
    else:
        pass
    return new_psi_vector, new_alpha
        
    
    
    

