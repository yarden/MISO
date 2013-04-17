##
## MISO scoring functions in Cython for MCMC sampler
##
import numpy as np
cimport numpy as np

np.import_array()
cimport cython

from libc.math cimport log

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_float_t

# cdef log_score_assignments(cnp.ndarray[DTYPE_t, ndim=1] isoform_nums,
#                            cnp.ndarray[double, ndim=1] psi_vector,
#                            cnp.ndarray[long, ndim=1] scaled_lens,
#                            int num_reads):
#     """
#     Single-end
#     ----------
#     Score an assignment of a set of reads given psi
#     and a gene (i.e. a set of isoforms).
#     """
#     cdef:
#        cnp.ndarray[double, ndim=1] psi_frag
#        cnp.ndarray[double, ndim=2] psi_frags
#     psi_frag = np.log(psi_vector) + np.log(scaled_lens)
#     psi_frag = psi_frag - scipy.misc.logsumexp(psi_frag)
#     psi_frags = np.tile(psi_frag, [num_reads, 1])
#     return psi_frags[np.arange(num_reads), isoform_nums]


##
## MCMC proposal functions
## 
# cdef propose_psi_vector(cnp.ndarray[DTYPE_float_t, ndim=1] alpha_vector):
#     """
#     Propose a new Psi vector.  Depends only on the alpha_vector
#     of parameters of the Dirichlet distribution from which the
#     current Psi vector was drawn.
#     """
#     cdef:
#         cnp.ndarray[DTYPE_float_t, ndim=1] proposed_psi_vector
#         cnp.ndarray[DTYPE_float_t, ndim=1] proposed_alpha_vector
#     proposed_psi_vector, proposed_alpha_vector = \
#         propose_norm_drift_psi_alpha(alpha_vector)
#     return (proposed_psi_vector, proposed_alpha_vector)




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def sum_log_score_reads(np.ndarray[DTYPE_t, ndim=2] reads,
                        np.ndarray[DTYPE_t, ndim=1] isoform_nums,
                        np.ndarray[DTYPE_t, ndim=1] num_parts_per_isoform,
                        np.ndarray[DTYPE_t, ndim=1] iso_lens,
                        int num_reads,
                        np.ndarray[double, ndim=1] log_num_reads_possible_per_iso):
    """
    Returns the sum of vectors of scores computed by 'log_score_reads'.
    """
    cdef np.ndarray[double, ndim=1] vect_log_scores
    cdef float sum_scores = 0
    cdef int curr_read = 0
    # Call log score reads to get vector of scores
    vect_log_scores = log_score_reads(reads,
                                      isoform_nums,
                                      num_parts_per_isoform,
                                      iso_lens,
                                      num_reads,
                                      log_num_reads_possible_per_iso)
    for curr_read in xrange(num_reads):
        sum_scores += vect_log_scores[curr_read]
    return sum_scores

cdef extern from "math.h":
    float INFINITY
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def log_score_reads(np.ndarray[DTYPE_t, ndim=2] reads,
                    np.ndarray[DTYPE_t, ndim=1] isoform_nums,
                    np.ndarray[DTYPE_t, ndim=1] num_parts_per_isoform,
                    np.ndarray[DTYPE_t, ndim=1] iso_lens,
                    int num_reads,
                    np.ndarray[double, ndim=1] log_num_reads_possible_per_iso):
    """
    Scores a set of reads given their isoform assignments.

    Parameters:
    -----------
    reads : 2-d array, representation of reads
    isoform_nums : 1-d np.array, assignment of each read to isoform. ith entry
                   is the assignment of ith read.
    num_parts_per_isoform : 1-d array, number of parts (exons) in each isoform
    iso_lens : 1-d array, lengths of each isoform
    num_reads : int, number of reads to process
    log_num_reads_possible_per_iso : 1-d array, the log'd number of reads
    possible in each isoform taking into account read length and the
    overhang length constraint.
    """
    cdef np.ndarray[double, ndim=1] log_prob_reads = np.empty(num_reads)
    # Read counter
    cdef int curr_read = 0
    # Isoform counter
    cdef int curr_iso_num = 0
    # Constant used in probability calculation
    cdef double log_one_val = log(1)
    #cdef double log_zero_prob_val = -2000
    #cdef double log_zero_prob_val = -1 * float('inf')
    cdef double log_zero_prob_val = -1 * INFINITY
    for curr_read in xrange(num_reads):
        # For each isoform assignment, score its probability
        # Get the current isoform's number (0,...,K-1 for K isoforms)
        curr_iso_num = isoform_nums[curr_read]
        # Check if the read is consistent with isoform
        # if it isn't, record 0 probability
        if reads[curr_read, curr_iso_num] == 0:
            # Read consistent with isoform
            log_prob_reads[curr_read] = log_zero_prob_val
        else:
            log_prob_reads[curr_read] = \
                log_one_val - log_num_reads_possible_per_iso[curr_iso_num]
    return log_prob_reads



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





