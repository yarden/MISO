##
## Proposal functions for Metropolis-Hastings sampler
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython

cimport stat_helpers
cimport matrix_utils
cimport math_utils
cimport sampling_utils
cimport array_utils

def propose_norm_drift_psi_alpha(double[:] alpha_vector,
                                 double[:, :] sigma_mat,
                                 double[:, :] L):
    """
    Propose a new alpha vector (parameters to Beta/Dirichlet distribution)
    and Psi vector (drawn from a Beta with those alpha vector)
    given a sigma matrix. The new alpha vector is drawn
    from a logistic normal distribution, i.e. by drawing
    a vector from a multivariate normal distribution and then
    taking an inverse logit transform (logit^-1()) to get
    values on the simplex.

    Args:
      alpha_vector: Vector of alpha parameters to Beta/Dirichlet
      distribution
      
      sigma_mat: Covariance matrix (2d array)

      L: Cholesky decomposition of covariance matrix (sigma_mat)

    Returns:
      psi vector (1d) and alpha vector (1d)

    if k is number of isoforms, psi vector is of length k
    and alpha vector is of length k-1
    """
    cdef int alpha_vect_len = alpha_vector.shape[0]
    # Psi vector without last entry (k-1 dimensions)
    cdef double[:] new_partial_psi = \
      array_utils.get_double_array(alpha_vect_len)
    cdef double partial_psi_sum = 0.0
    # Complete Psi vector
    cdef double[:] new_psi_vector = \
      array_utils.get_double_array(alpha_vect_len + 1)
    cdef double[:] new_alpha = \
      array_utils.get_double_array(alpha_vect_len)
    cdef int n = 0
    # First entry in covariance matrix
    #cdef double covar_first_entry = sigma_mat[0, 0]
    # Sample new alpha vector
    new_alpha = \
      sampling_utils.sample_multivar_normal(alpha_vector, L, alpha_vect_len)
    # Do the inverse logit transform to get a set of Psi vectors of
    # dimension k - 1, where k is number of isoforms (i.e. missing the
    # last entry)
    new_partial_psi = sampling_utils.logit_inv(new_alpha, alpha_vect_len)
    partial_psi_sum = array_utils.sum_array(new_partial_psi, alpha_vect_len)
    # Set entries of the new psi vector to be that of
    # the Psi values we get for the first k-1 entries
    # through the inverse logit transform
    for n in xrange(alpha_vect_len):
        new_psi_vector[n] = new_partial_psi[n]
    # Set the last entry of vector to be 1 - sum(partial_psi)
    new_psi_vector[alpha_vect_len + 1] = 1 - partial_psi_sum
    return new_psi_vector, new_alpha

    
    
    

