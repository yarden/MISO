##
## MISO scoring functions in Cython for MCMC sampler
##
import numpy as np
cimport numpy as np

cnp.import_array()
cimport cython

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_float_t

cdef log_score_reads(np.ndarray[DTYPE_t, ndim=2] reads,
                     np.ndarray[DTYPE_t, ndim=1] isoform_nums,
                     np.ndarray[DTYPE_t, ndim=1] num_parts_per_isoform,
                     np.ndarray[DTYPE_t, ndim=1] iso_lens,
                     int read_len,
                     int overhang_len,
                     int num_reads):
    """
    Single-end
    ----------
    Score a set of reads given their isoform assignments.
    """
    cdef:
       cnp.ndarray[double, ndim=1] log_prob_reads
       cnp.ndarray[long, ndim=1] overhang_excluded
       cnp.ndarray[long, ndim=1] zero_prob_indx
       cnp.ndarray[long, ndim=1] num_reads_possible
    # The probability of a read being assigned to an isoform that
    # could not have generated it (i.e. where the read is not a
    # substring of the isoform) is zero.  Check for consistency
    overhang_excluded = \
        2*(overhang_len - 1)*(num_parts_per_isoform[isoform_nums] - 1)
    # The number of reads possible is the number of ways a 36 nt long read can be
    # aligned onto the length of the current isoforms, minus the number of positions
    # that are ruled out due to overhang constraints.
    num_reads_possible = \
        (iso_lens[isoform_nums] - read_len + 1) - overhang_excluded
    log_prob_reads = np.log(1) - np.log(num_reads_possible)
    zero_prob_indx = np.nonzero(reads[np.arange(num_reads), isoform_nums] == 0)[0]
    # Assign probability 0 to reads inconsistent with assignment
    log_prob_reads[zero_prob_indx] = -1 * np.inf
    return log_prob_reads


cdef log_score_assignments(cnp.ndarray[DTYPE_t, ndim=1] isoform_nums,
                           cnp.ndarray[double, ndim=1] psi_vector,
                           cnp.ndarray[long, ndim=1] scaled_lens,
                           int num_reads):
    """
    Single-end
    ----------
    Score an assignment of a set of reads given psi
    and a gene (i.e. a set of isoforms).
    """
    cdef:
       cnp.ndarray[double, ndim=1] psi_frag
       cnp.ndarray[double, ndim=2] psi_frags
    psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    psi_frag = psi_frag - scipy.misc.logsumexp(psi_frag)
    psi_frags = np.tile(psi_frag, [num_reads, 1])
    return psi_frags[np.arange(num_reads), isoform_nums]


##
## MCMC proposal functions
## 
cdef propose_psi_vector(cnp.ndarray[DTYPE_float_t, ndim=1] alpha_vector):
    """
    Propose a new Psi vector.  Depends only on the alpha_vector
    of parameters of the Dirichlet distribution from which the
    current Psi vector was drawn.
    """
    cdef:
        cnp.ndarray[DTYPE_float_t, ndim=1] proposed_psi_vector
        cnp.ndarray[DTYPE_float_t, ndim=1] proposed_alpha_vector
    proposed_psi_vector, proposed_alpha_vector = \
        propose_norm_drift_psi_alpha(alpha_vector)
    return (proposed_psi_vector, proposed_alpha_vector)


#     def propose_norm_drift_psi_alpha(self, alpha_vector):
#         if len(alpha_vector) == 1:
#             alpha_vector = alpha_vector[0]
# #            print "proposing from normal with mean: ", alpha_vector, " exp: ", exp(alpha_vector)
#             alpha_next = [normal(alpha_vector, self.params['sigma_proposal'])]
# #            print "got alpha_next: ", alpha_next, " exp: ", exp(alpha_next)
#             new_psi = logit_inv([alpha_next[0]])[0]
#             new_psi_vector = [new_psi, 1-new_psi]
#         else:
#             alpha_next = multivariate_normal(alpha_vector, self.params['sigma_proposal'])
#             new_psi = logit_inv(alpha_next)
#             new_psi_vector = concatenate((new_psi, array([1-sum(new_psi)])))
#         return (new_psi_vector, alpha_next)





