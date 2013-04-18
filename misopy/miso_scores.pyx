##
## MISO scoring functions in Cython for MCMC sampler
##
import numpy as np
cimport numpy as np

np.import_array()
cimport cython

from libc.math cimport log, exp

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_float_t


##
## Statistics helper functions (TODO: remove to separate module)
##
cdef my_logsumexp(np.ndarray[DTYPE_float_t, ndim=1] log_vector,
                  int vector_len):
    """
    Log sum exp.

    Parameters:
    -----------

    log_vector : array of floats corresponding to log values.
    vector_len : int, length of vector.

    Returns:
    --------

    Result of log(sum(exp(log_vector)))
    """
    cdef DTYPE_float_t curr_exp_value = 0.0
    cdef DTYPE_float_t sum_of_exps = 0.0
    cdef DTYPE_float_t log_sum_of_exps = 0.0
    cdef int curr_elt = 0
    for curr_elt in xrange(vector_len):
        curr_exp_value = exp(log_vector[curr_elt])
        sum_of_exps += curr_exp_value
    # Now take log of the sum of exp values
    log_sum_of_exps = log(sum_of_exps)
    return log_sum_of_exps


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


##
## Sampler scoring functions
##
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def log_score_assignments(np.ndarray[DTYPE_t, ndim=1] isoform_nums,
                          np.ndarray[DTYPE_float_t, ndim=1] log_psi_frag_vector,
                          int num_reads):
    """
    Score an assignment of a set of reads given psi
    and a gene (i.e. a set of isoforms).
    """
    cdef np.ndarray[DTYPE_float_t, ndim=1] log_scores = np.empty(num_reads)
    cdef DTYPE_float_t curr_log_psi_frag = 0.0
    cdef int curr_read = 0
    cdef int curr_iso_num = 0
    for curr_read in xrange(num_reads):
        curr_iso_num = isoform_nums[curr_read]
        curr_log_psi_frag = log_psi_frag_vector[curr_iso_num]
        log_scores[curr_read] = curr_log_psi_frag 
    return log_scores


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def sum_log_score_assignments(np.ndarray[DTYPE_t, ndim=1] isoform_nums,
                              np.ndarray[DTYPE_float_t, ndim=1] log_psi_frag_vector,
                              int num_reads):
    """
    Score an assignment of a set of reads given psi
    and a gene (i.e. a set of isoforms).
    """
    cdef np.ndarray[DTYPE_float_t, ndim=1] assignment_scores = \
        np.empty(num_reads)
    cdef DTYPE_float_t sum_log_scores = 0.0
    cdef DTYPE_float_t curr_assignment_score = 0.0
    cdef int curr_read = 0
    # Get log score of assignments
    assignment_scores = log_score_assignments(isoform_nums,
                                              log_psi_frag_vector,
                                              num_reads)
    for curr_read in xrange(num_reads):
        curr_assignment_score = assignment_scores[curr_read]
        sum_log_scores += curr_assignment_score
    return sum_log_scores


def compute_log_psi_frag(np.ndarray[DTYPE_float_t, ndim=1] psi_vector,
                         np.ndarray[DTYPE_t, ndim=1] scaled_lens,
                         int num_isoforms):
    """
    Compute log Psi frag from Psi vector.

    FOR SINGLE-END right now.
    """
    # Log psi frag vector computed from psi vector
    cdef np.ndarray[DTYPE_float_t, ndim=1] log_psi_frag = \
        np.empty(num_isoforms)
    # Isoform counter
    cdef int curr_isoform = 0
    for curr_isoform in xrange(num_isoforms):
        log_psi_frag[curr_isoform] = \
            log(psi_vector[curr_isoform]) + log(scaled_lens[curr_isoform])
    # Normalize scaled Psi values to sum to 1
    log_psi_frag = \
        log_psi_frag - my_logsumexp(log_psi_frag, num_isoforms)
    return log_psi_frag
    #curr_log_psi_frag = \
    #    np.log(curr_psi_vector) + np.log(self.scaled_lens_single_end)
    #curr_log_psi_frag = \
    #    curr_log_psi_frag - scipy.misc.logsumexp(curr_log_psi_frag)


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
    ##
    ## TODO: move 'log_one_val' and 'zero_log_prob_val' to be
    ## global constants
    ##
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





