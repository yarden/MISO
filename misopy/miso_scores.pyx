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

# /* We only handle a special case here, where sigma is a diagonal
#    matrix with identical elements. In this case it is easy to 
#    invert it, or calculate its determinant. */

# int splicing_mvplogisnorm(const splicing_vector_t *theta, 
# 			  const splicing_vector_t *mu, 
# 			  double sigma, int len, double *score) {
  
#   double covarConst = pow(2 * M_PI * sigma, -0.5 * len);
#   double ltheta=1.0, prodTheta=1.0, expPart=0.0, pdfVal;
#   int i;
  
#   for (i=0; i<len; i++) {
#     double at=VECTOR(*theta)[i];
#     ltheta -= at;
#     prodTheta *= at;
#   }
#   prodTheta = 1.0 / prodTheta / ltheta;
  
#   for (i=0; i<len; i++) {
#     double tmp=log(VECTOR(*theta)[i] / ltheta) - VECTOR(*mu)[i];
#     expPart += (-0.5) * tmp * tmp / sigma;
#   }
  
#   pdfVal = covarConst * prodTheta * exp(expPart);

#   *score = log(pdfVal);
  
#   return 0;
# }



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def my_cumsum(np.ndarray[double, ndim=1] input_array):
    """
    Return cumulative sum of array.
    """
    # Cumulative sum at every position
    cdef np.ndarray[double, ndim=1] cumsum_array = np.empty_like(input_array)
    cdef int curr_elt = 0
    cdef int num_elts = input_array.shape[0]
    # Current cumulative sum: starts at first element
    cdef double curr_cumsum = 0.0
    for curr_elt in xrange(num_elts):
        cumsum_array[curr_elt] = (input_array[curr_elt] + curr_cumsum)
        curr_cumsum = cumsum_array[curr_elt]
    return cumsum_array


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef DTYPE_float_t my_logsumexp(np.ndarray[DTYPE_float_t, ndim=1] log_vector,
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
## MCMC sampler logic functions
##
# def compute_metropolis_ratio(self, reads, assignments, proposed_psi_vector,
#                              proposed_alpha_vector,
#                              curr_psi_vector, curr_alpha_vector,
#                              hyperparameters,
#                              full_metropolis):
#     """
#     Compute the Metropolis-Hastings ratio:

#         P(psi_next)Q(psi; psi_next)
#         ---------------------------
#         P(psi)Q(psi_next; psi)

#     Parameters:
#     -----------
#     reads : array, reads to be processed
#     assignments : array, assignments of reads to isoforms
#     proposed_psi_vector : array, proposed Psi vector
#     proposed_alpha_vector : array, proposed alpha vector (parameter to Dirichlet)
#     from which the Psi value vector was drawn
#     curr_psi_vector : array, current Psi vector
#     curr_alpha_vector : array, current alpha vector (parameter to Dirichlet)
#     hyperparameters : array, hyperparameters for scoring
#     full_metropolis : int, if 1 then compute full MH ratio, otherwise not
#     """
#     # Compute acceptance ratio: the joint score for proposed Psi divided
#     # by joint score given current Psi
#     # P(Psi', ...)
#     proposed_joint_score = \
#         self.log_score_joint(reads, assignments, proposed_psi_vector,
#                              gene, hyperparameters)
#     # P(Psi, ...)
#     curr_joint_score = \
#         self.log_score_joint(reads, assignments, curr_psi_vector,
#                              gene, hyperparameters)
#     if curr_joint_score == -inf:
#         self.miso_logger.error("Joint score of current state is negative infinity!")
#         self.miso_logger.error("  - assignments: " + str(assignments))
#         self.miso_logger.error("  - psi vector: " + str(curr_psi_vector))
#         self.miso_logger.error("  - reads: " + str(reads))
#         raise Exception, "curr_joint_score is negative."
#     # Q(x; x'), the probability of proposing to move back to current state from
#     # proposed state x'
#     mh_ratio = None
#     proposal_to_curr_score = \
#         log_score_psi_vector_transition(curr_psi_vector, proposed_alpha_vector)
#     # Q(x'; x), the probability of proposing to move to the proposed state x' from
#     # the current state
#     curr_to_proposal_score = \
#         log_score_psi_vector_transition(proposed_psi_vector, curr_alpha_vector)
#     # Computing full Metropolis-Hastings ratio
#     if full_metropolis == 0:
#         # Not full MH; just ratio of proposed to current joint score
#         mh_ratio = (proposed_joint_score - curr_joint_score)
#     else:
#         # Full MH ratio
#         mh_ratio = (proposed_joint_score + proposal_to_curr_score) - \
#                    (curr_joint_score + curr_to_proposal_score)
#     if curr_to_proposal_score == (-1 * INFINITY):
#         self.miso_logger.error("curr to proposal is -inf")
#         raise Exception, "curr to proposal is -Inf"
#     if proposed_joint_score == (-1 * INFINITY):
#         self.miso_logger.debug("Proposing to move to impossible state!")	    
#         raise Exception, "Proposing to move to impossible state!"
#     if abs(mh_ratio) == Inf:
#         self.miso_logger.debug("MH ratio is Inf!")
#         raise Exception, "MH ratio is Inf!"
#     return (exp(mh_ratio), curr_joint_score, proposed_joint_score)    

    

##
## Sampler scoring functions
##
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def dirichlet_lnpdf(np.ndarray[double, ndim=1] alpha,
                    np.ndarray[double, ndim=1] vector):
    """
    Wrapper for dirichlet log pdf scoring function.
    """
    cdef int D = vector.size
    return dirichlet_log_pdf_raw(D,
                                 &alpha[0], alpha.strides[0],
                                 &vector[0], vector.strides[0])


@cython.infer_types(True)
cdef double dirichlet_log_pdf_raw(
    int D,
    double* alpha, int alpha_stride,
    double* vector, int vector_stride,
    ):
    """Compute the log of the Dirichlet PDF evaluated at one vector."""

    cdef void* alpha_p = alpha
    cdef void* vector_p = vector

    # first term
    term_a = 0.0

    for d in xrange(D):
        term_a += (<double*>(alpha_p + alpha_stride * d))[0]

    term_a = libc.math.lgamma(term_a)

    # second term
    term_b = 0.0

    for d in xrange(D):
        term_b += libc.math.lgamma((<double*>(alpha_p + alpha_stride * d))[0])

    # third term
    cdef double alpha_d
    cdef double vector_d

    term_c = 0.0

    for d in xrange(D):
        alpha_d = (<double*>(alpha_p + alpha_stride * d))[0]
        vector_d = (<double*>(vector_p + vector_stride * d))[0]

        term_c += (alpha_d - 1.0) * libc.math.log(vector_d)

    # ...
    return term_a - term_b + term_c

# def log_score_joint_single_end(np.ndarray[DTYPE_t, ndim=1] reads,
#                                np.ndarray[DTYPE_t, ndim=1] assignments,
#                                np.ndarray[DTYPE_float_t, ndim=1] psi_vector,
#                                DTYPE_float_t log_score_psi_vector,
#                                int num_reads):
#     """
#     Return a log score for the joint distribution for single-end reads.

#     Parameters:
#     -----------

#     reads : array, array of reads
#     assignments : assignments of reads to isoforms
#     psi_vector : array, Psi vector in current state
#     log_score_psi_vector : float, the log score of Psi vector (computed
#     from Python)
#     num_reads : int, number of reads
#     """
#     # Get the logged Psi frag vector
#     # ....
#     # Score the read
#     if not self.paired_end:
# #            log_reads_prob = \
# #                sum(self.log_score_reads(reads, assignments, gene))
# #	    log_reads_prob = \
# #                sum(miso_scores.log_score_reads(reads,
# #                                                assignments,
# #                                                self.num_parts_per_isoform,
# #                                                self.iso_lens,
# #                                                self.num_reads,
# #                                                self.log_num_reads_possible_per_iso))
#         log_reads_prob = \
#             miso_scores.sum_log_score_reads(reads,
#                                             assignments,
#                                             self.num_parts_per_isoform,
#                                             self.iso_lens,
#                                             self.num_reads,
#                                             self.log_num_reads_possible_per_iso)
#         log_psi_frag = \
#             miso_scores.compute_log_psi_frag(psi_vector,
#                                              self.scaled_lens_single_end,
#                                              self.num_isoforms)
#     else:
#         raise Exception, "Paired-end not implemented!"
#         log_reads_prob = \
#             sum(self.log_score_paired_end_reads(reads, assignments, gene))
#         log_psi_frag = None
#     if not self.paired_end:
#         log_assignments_prob = \
#             miso_scores.sum_log_score_assignments(assignments,
#                                                   log_psi_frag,
#                                                   self.num_reads)
#     else:
#         raise Exception, "Paired-end not implemented!"
#         log_assignments_prob = \
#             sum(self.log_score_paired_end_assignment(assignments,
#                                                      psi_vector,
#                                                      gene))
#     # Score the Psi vector: keep this in Python for now
#     log_psi_prob = self.log_score_psi_vector(psi_vector, hyperparameters)
#     log_joint_score = log_reads_prob + log_assignments_prob + log_psi_prob
#     return log_joint_score


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def log_score_psi_vector(np.ndarray[DTYPE_float_t, ndim=1] psi_vector,
                         np.ndarray[DTYPE_float_t, ndim=1] hyperparameters):
    return dirichlet_lnpdf(hyperparameters, psi_vector)


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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
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





