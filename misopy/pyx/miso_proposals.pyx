##
## Proposal functions and related functions for Metropolis-Hastings sampler
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython

from libc.math cimport log
from libc.math cimport exp
from libc.stdlib cimport rand

cdef float MY_MAX_INT = float(10000)
# Define infinity
cdef extern from "math.h":
    float INFINITY

NEG_INFINITY = (-1 * INFINITY)

cimport stat_helpers
cimport matrix_utils
cimport math_utils
cimport sampling_utils
cimport array_utils
cimport miso_scores_single as scores_single

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
    new_partial_psi = math_utils.logit_inv(new_alpha, alpha_vect_len)
    partial_psi_sum = array_utils.sum_array(new_partial_psi, alpha_vect_len)
    # Set entries of the new psi vector to be that of
    # the Psi values we get for the first k-1 entries
    # through the inverse logit transform
    for n in xrange(alpha_vect_len):
        new_psi_vector[n] = new_partial_psi[n]
    # Set the last entry of vector to be 1 - sum(partial_psi)
    new_psi_vector[alpha_vect_len] = 1 - partial_psi_sum
    return new_psi_vector, new_alpha


cdef double compute_mh_ratio(int[:, :] reads,
                             int[:] frag_lens,
                             int[:] curr_assignments,
                             int[:] new_assignments,
                             double[:] curr_psi_vector,
                             double[:] new_psi_vector,
                             double[:] curr_log_psi_frag,
                             double[:] new_log_psi_frag,
                             double[:] curr_assignment_scores,
                             double[:] new_assignment_scores,
                             int[:] num_parts_per_isoform,
                             int[:] iso_lens,
                             double[:] log_num_reads_possible_per_iso,
                             int read_len,
                             int num_reads,
                             double[:] hyperparameters,
                             int overhang_len,
                             int paired_end,
                             int full_metropolis):
    """
    Compute the Metropolis-Hastings ratios.

        P(psi_next)Q(psi; psi_next)
        ---------------------------
        P(psi)Q(psi_next; psi)

    Args:
      ...

    Kwargs:
      ...

    Returns:
      ...
    """
    ###
    ### TODO: Handle paired_end parameter here!
    ###
    # Compute acceptance ratio: the joint score for proposed Psi divided
    # by joint score given current Psi
    # P(Psi', ...)
    proposed_joint_score = \
      scores_single.log_score_joint_single_end(reads,
                                               curr_assignments,
                                               curr_psi_vector,
                                               curr_log_psi_frag,
                                               curr_assignment_scores,
                                               num_parts_per_isoform,
                                               iso_lens,
                                               log_num_reads_possible_per_iso,
                                               read_len,
                                               num_reads,
                                               hyperparameters,
                                               overhang_len)      
    # P(Psi, ...)
    curr_joint_score = \
        scores_single.log_score_joint(reads,
                                      new_assignments,
                                      new_psi_vector,
                                      new_log_psi_frag,
                                      new_assignment_scores,
                                      num_parts_per_isoform,
                                      iso_lens,
                                      log_num_reads_possible_per_iso,
                                      read_len,
                                      num_reads,
                                      hyperparameters,
                                      overhang_len)      
    # if curr_joint_score == -inf:
    #     self.miso_logger.error("Joint score of current state is negative infinity!")
    #     self.miso_logger.error("  - assignments: " + str(assignments))
    #     self.miso_logger.error("  - psi vector: " + str(curr_psi_vector))
    #     self.miso_logger.error("  - reads: " + str(reads))
    #     raise Exception, "curr_joint_score is negative."
    # Q(x; x'), the probability of proposing to move back to current state from
    # proposed state x'
    mh_ratio = None
    # Score logistic normal here
    proposal_to_curr_score = 0.0
    # Q(x'; x), the probability of proposing to move to the proposed 
    # state x' from the current state
    curr_to_proposal_score = 0.0
    # Computing full Metropolis-Hastings ratio
    if not full_metropolis:
        mh_ratio = (proposed_joint_score - curr_joint_score)
    else:
        mh_ratio = (proposed_joint_score + proposal_to_curr_score) - \
                   (curr_joint_score + curr_to_proposal_score)
    if curr_to_proposal_score == NEG_INFINITY:
        raise Exception, "curr to proposal is -Inf"
    if proposed_joint_score == NEG_INFINITY:
        raise Exception, "Proposing to move to impossible state!"
    if abs(mh_ratio) == INFINITY:
        raise Exception, "MH ratio is Inf!"
    return (exp(mh_ratio), curr_joint_score, proposed_joint_score)
