##
## Proposal functions for Metropolis-Hastings sampler
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython

cimport stat_helpers
cimport matrix_utils
cimport math_utils
cimport array_utils

cpdef propose_norm_drift_psi_alpha(double[:] alpha_vector,
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
    cdef double[:] new_psi_vector = \
      array_utils.get_double_array(alpha_vect_len + 1)
    cdef double[:] new_alpha = \
      array_utils.get_double_array(alpha_vect_len)
    if alpha_vect_len == 1:
        pass
    else:
        pass
    return new_psi_vector, new_alpha


    # def propose_norm_drift_psi_alpha(self, alpha_vector):
    #     if len(alpha_vector) == 1:
    #         alpha_vector = alpha_vector[0]
    #         alpha_next = [normal(alpha_vector, self.params['sigma_proposal'])]
    #         new_psi = logit_inv([alpha_next[0]])[0]
    #         new_psi_vector = [new_psi, 1-new_psi]
    #     else:
    #         alpha_next = multivariate_normal(alpha_vector, self.params['sigma_proposal'])
    #         new_psi = logit_inv(alpha_next)
    #         new_psi_vector = concatenate((new_psi, array([1-sum(new_psi)])))
    #     return (new_psi_vector, alpha_next)

        
    
    
    

