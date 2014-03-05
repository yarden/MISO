##
## MISO scoring functions in Cython for MCMC sampler
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython
from cython.view cimport array as cvarray
from cpython.array cimport array, clone

cimport stat_helpers
cimport matrix_utils
cimport array_utils
cimport math_utils

from libc.math cimport log
from libc.math cimport exp
from libc.stdlib cimport rand

cdef float MY_MAX_INT = float(10000)
# Define infinity
cdef extern from "math.h":
    float INFINITY

NEG_INFINITY = (-1 * INFINITY)

import misopy


##
## Sampler scoring functions
##
cpdef double log_score_joint_single_end(int[:, :] reads,
                                        int[:] isoform_nums,
                                        double[:] psi_vector,
                                        double[:] log_psi_frag,
                                        double[:] assignment_scores,
                                        int[:] num_parts_per_isoform,
                                        int[:] iso_lens,
                                        double[:] log_num_reads_possible_per_iso,
                                        int read_len,
                                        int num_reads,
                                        double[:] hyperparameters,
                                        int overhang_len):
    """
    Return a log score for the joint distribution for single-end reads.

    The joint score requires scoring:

      P(reads|assignments) x P(assignments|Psi) x P(Psi)

    This computes the sum of log scores of these:

      logP(reads|assignments) + logP(assignments|Psi) + logP(psi)
    
    Args:
      reads: the read alignments (2d, int array)
      isoform_nums: the assignments of reads to isoforms
      psi_vector: the psi vector to consider 
      num_reads: number of reads

    Returns:
      log probability score
    """
    cdef double log_reads_prob
    cdef double log_assignments_prob
    cdef double log_psi_prob
    cdef double log_joint_score
    # log P(reads|assignments)
    log_reads_prob = \
      sum_log_score_reads(reads, isoform_nums, num_parts_per_isoform,
                          iso_lens, log_num_reads_possible_per_iso,
                          num_reads, read_len, overhang_len)
    # log P(assignments|Psi)
    log_assignments_prob = \
      sum_log_score_assignments(isoform_nums,
                                log_psi_frag,
                                num_reads,
                                assignment_scores)
    # log P(Psi)
    log_psi_prob = log_score_psi_vector(psi_vector, hyperparameters)
    # Joint score
    log_joint_score = log_reads_prob + log_assignments_prob + log_psi_prob
    return log_joint_score


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef double log_score_psi_vector(double[:] psi_vector,
                                  double[:] hyperparameters):
    return stat_helpers.dirichlet_lnpdf(hyperparameters, psi_vector)


#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)
cpdef double[:] \
  log_score_assignments(int[:] isoform_nums,
                        double[:] log_psi_frag_vector,
                        int num_reads,
                        double[:] log_scores):
    """
    Score an assignment of a set of reads given psi
    and a gene (i.e. a set of isoforms).
    """
    cdef double curr_log_psi_frag = 0.0
    cdef int curr_read = 0
    cdef int curr_iso_num = 0
    for curr_read in xrange(num_reads):
        curr_iso_num = isoform_nums[curr_read]
        # The score of an assignment to isoform i is the ith entry
        # is simply the Psi_Frag vector
        curr_log_psi_frag = log_psi_frag_vector[curr_iso_num]
        log_scores[curr_read] = curr_log_psi_frag 
    return log_scores


#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)
cpdef double \
  sum_log_score_assignments(int[:] isoform_nums,
                            double[:] log_psi_frag_vector,
                            int num_reads,
                            double[:] assignment_scores):
    """
    Score an assignment of a set of reads given psi
    and a gene (i.e. a set of isoforms).
    """
    cdef double sum_log_scores = 0.0
    cdef double curr_assignment_score = 0.0
    cdef int curr_read = 0
    # Get log score of assignments
    assignment_scores = log_score_assignments(isoform_nums,
                                              log_psi_frag_vector,
                                              num_reads,
                                              assignment_scores)
    for curr_read in xrange(num_reads):
        curr_assignment_score = assignment_scores[curr_read]
        sum_log_scores += curr_assignment_score
    return sum_log_scores


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef double[:] compute_log_psi_frag(double[:] psi_vector,
                                     double[:] scaled_lens,
                                     int num_isoforms):
    """
    Compute log Psi frag from Psi vector.

    FOR SINGLE-END right now.
    """
    # Log psi frag vector computed from psi vector
    cdef log_psi_frag = array_utils.get_double_array(num_isoforms)
    # Isoform counter
    cdef int curr_isoform = 0
    for curr_isoform in xrange(num_isoforms):
        log_psi_frag[curr_isoform] = \
            log(psi_vector[curr_isoform]) + log(scaled_lens[curr_isoform])
    # Normalize scaled Psi values to sum to 1
    log_psi_frag = \
        log_psi_frag - math_utils.my_logsumexp(log_psi_frag, num_isoforms)
    return log_psi_frag


#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)
cpdef double sum_log_score_reads(int[:, :] reads,
                                 int[:] isoform_nums,
                                 int[:] num_parts_per_isoform,
                                 int[:] iso_lens,
                                 double[:] log_num_reads_possible_per_iso,
                                 int num_reads,
                                 int read_len,
                                 int overhang_len):
    """
    Returns the sum of vectors of scores computed by 'log_score_reads'.
    """
    cdef double[:] vect_log_scores = \
      array_utils.get_double_array(num_reads)
    cdef double sum_scores = 0.0
    cdef int curr_read = 0
    # Call log score reads to get vector of scores
    vect_log_scores = log_score_reads(reads,
                                      isoform_nums,
                                      num_parts_per_isoform,
                                      iso_lens,
                                      log_num_reads_possible_per_iso,
                                      num_reads,
                                      read_len,
                                      overhang_len,
                                      vect_log_scores)
    for curr_read in xrange(num_reads):
        # If a score for any of the reads is -inf, then
        # the sum of the scores is -inf, so no point computing
        # the rest
        if vect_log_scores[curr_read] == (-1 * INFINITY):
            sum_scores = (-1 * INFINITY)
            break
        # Sum the scores
        sum_scores += vect_log_scores[curr_read]
    return sum_scores

    
#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)
cpdef double[:] log_score_reads(int[:, :] reads,
                                int[:] isoform_nums,
                                int[:] num_parts_per_iso,
                                int[:] iso_lens,
                                double[:] log_num_reads_possible_per_iso,
                                int num_reads,
                                int read_len,
                                int overhang_len,
                                double[:] log_prob_reads):
    """
    Scores a set of reads given their isoform assignments.

    Parameters:
    -----------
    reads : 2-d array, representation of reads
    isoform_nums : 1-d np.array, assignment of each read to isoform. ith entry
                   is the assignment of ith read.
    num_parts_per_iso : 1-d array, number of parts (exons) in each isoform
    iso_lens : 1-d array, lengths of each isoform
    num_reads : int, number of reads to process
    log_num_reads_possible_per_iso : 1-d array, the log'd number of reads
    possible in each isoform taking into account read length and the
    overhang length constraint.
    """
    # Read counter
    cdef int curr_read = 0
    # Isoform counter
    cdef int curr_iso_num = 0
    for curr_read in xrange(num_reads):
        # For each isoform assignment, score its probability
        # Get the current isoform's number (0,...,K-1 for K isoforms)
        curr_iso_num = isoform_nums[curr_read]
        # Check if the read is consistent with isoform
        # if it isn't, record 0 probability
        if reads[curr_read, curr_iso_num] == 0:
            # Read consistent with isoform
            log_prob_reads[curr_read] = NEG_INFINITY
        else:
            # Uniform scoring: probability of read is:
            #    1/(number of possible reads from isoform)
            #  = 1/(iso_len - read_len + 1)
            #  = log[1/(iso_len - read_len + 1)]
            #  = log(1) - log(iso_len - read_len + 1)
            log_prob_reads[curr_read] = \
                log(1) - log_num_reads_possible_per_iso[curr_iso_num]
    return log_prob_reads


##
## Sampling functions
##
# Sample reassignments
#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)
#int[:, :] reads,
cpdef int[:] \
  sample_reassignments(int[:, :] reads,
                       double[:] psi_vector,
                       double[:] log_psi_frag_vector,
                       double[:] log_num_reads_possible_per_iso,
                       int[:] scaled_lens,
                       int[:] iso_lens,
                       int[:] num_parts_per_iso,
                       int[:] iso_nums,
                       int num_reads,
                       int read_len,
                       int overhang_len,
                       int[:] new_assignments):
    """
    Sample a reassignments of reads to isoforms.
    Note that this does not depend on the read's current assignment since
    we're already considering the possibility of 'reassigning' the read to
    its current assignment in the probability calculations.
    """
    cdef int num_isoforms = psi_vector.shape[0]
    # Initial probabilities of assignments given Psi vector
    cdef double[:] log_assignment_probs = \
      array_utils.get_double_array(num_reads)
    # Initial probabilties of reads given assignments
    cdef double[:] log_read_probs = \
      array_utils.get_double_array(num_reads)
    # Probabilities of reassigning current read to each of the isoforms
    cdef double[:] reassignment_probs = \
      array_utils.get_double_array(num_isoforms)
    # Normalized reassignment probabilities
    cdef double[:] norm_reassignment_probs = \
      array_utils.get_double_array(num_isoforms)
    # The candidate assignment of current read
    cdef int[:] cand_assignment = \
      array_utils.get_int_array(1)
    # Log probability of candidate assignment
    cdef double cand_assignment_log_prob
    # Array containing a single-read, the current read to be
    # considered
    ##### TODO: Optimize
    cdef int[:, :] curr_read_array = \
      cvarray(shape=(2,2), itemsize=sizeof(int), format="i")
    # Temporary vector to score assignment log probs 
    cdef double[:] temp_assignment_log_prob = \
      array_utils.get_double_array(1)
    cdef double[:] temp_reads_log_prob = \
      array_utils.get_double_array(1)
    cdef int[:] temp_curr_read_assignment = \
      array_utils.get_int_array(1)
    # The sampled current read assignment
    cdef int curr_read_assignment = 0
    # Current isoform counter
    cdef int curr_isoform = 0
    # Current read counter
    cdef int curr_read = 0
    # First calculate the scores of all the current reads
    # given their isoform assignment
    log_read_probs = log_score_reads(reads,
                                     iso_nums,
                                     num_parts_per_iso,
                                     iso_lens,
                                     log_num_reads_possible_per_iso,
                                     num_reads,
                                     read_len,
                                     overhang_len,
                                     log_read_probs)
    # Also calculate the scores of all current assignments
    log_assignment_probs = log_score_assignments(iso_nums,
                                                 log_psi_frag_vector,
                                                 num_reads,
                                                 log_assignment_probs)
    # For each read, compute the probability of reassigning it to
    # each of the isoforms and sample a new reassignment
    for curr_read in xrange(num_reads):
        # Copy the current assignment of read probability
        old_read_prob = log_read_probs[curr_read]
        # Copy the current assignment of read to isoform
        old_assignment = iso_nums[curr_read]
        # Copy the current probability of assignment
        old_assignment_log_prob = log_assignment_probs[curr_read]
        # Compute the probability of assigning the current read
        # to each isoform
        for curr_isoform in xrange(num_isoforms):
            # Now consider reassigning this read to the current isoform
            iso_nums[curr_read] = curr_isoform
            # Score this assignment of reads to isoforms
            cand_assignment[0] = curr_isoform
            # Compute probability of candidate assignment. This depends
            # only on the Psi Frag vectors:
            #    P(I(j,k) | PsiFrag)
            cand_assignment_log_prob = \
              log_score_assignments(cand_assignment,
                                    log_psi_frag_vector,
                                    1,
                                    temp_assignment_log_prob)[0]
            # Compute the probability of the reads given the current
            # assignment. The only term that changes is the probability
            # of the current read, since we changed its assignment
            #    P(R | I)
            curr_read_array = reads[curr_read:curr_read+1]
            # Score probability of current read given candidate
            # assignment
            cand_read_log_prob = \
              log_score_reads(curr_read_array,
                              cand_assignment,
                              num_parts_per_iso,
                              iso_lens,
                              log_num_reads_possible_per_iso,
                              1,
                              read_len,
                              overhang_len,
                              temp_reads_log_prob)[0]
            # Set the new probability of current read
            log_read_probs[curr_read] = cand_read_log_prob
            # Set the new probability of read's assignment
            log_assignment_probs[curr_read] = cand_assignment_log_prob
            # Now compute reassignment probabilities
            # Probability of assigning the read to this isoform
            # is the sum: log(P(reads | assignment)) + log(P(assignment | Psi))
            # P(log(reads | assignment) is the sum of the read scores
            # vector with the current read's reassignment
            reassignment_probs[curr_isoform] = \
              (array_utils.sum_array(log_assignment_probs, num_reads) + \
               array_utils.sum_array(log_read_probs, num_reads))
            # Copy the old assignment of the read to the isoform
            iso_nums[curr_read] = old_assignment
            # Copy the old assignment probability
            log_assignment_probs[curr_read] = old_assignment_log_prob
            # Copy the old probability of read given the isoform
            log_read_probs[curr_read] = old_read_prob
        # Normalize the reassignment probabilities
        norm_reassignment_probs = \
          norm_log_probs(reassignment_probs, num_isoforms)
        curr_read_assignment = \
          sample_from_multinomial(norm_reassignment_probs,
                                  1,
                                  temp_curr_read_assignment)[0]
        new_assignments[curr_read] = curr_read_assignment
    return new_assignments


cdef double[:] \
     norm_log_probs(double[:] log_probs,
                    int vector_len):
    """
    Normalize log probabilities to sum to 1. Returns an UNLOGGED
    probability!
    """
    # Normalizing factor
    cdef double norm_factor = math_utils.my_logsumexp(log_probs, vector_len)
    # Normalized probabilities (UNLOGGED)
    cdef double[:] norm_probs = array_utils.get_double_array(vector_len)
    for curr_entry in xrange(vector_len):
        # If the current log probability is -inf, then
        # record it as 0 (unlogged probability)
        if log_probs[curr_entry] == NEG_INFINITY:
            norm_prob = 0.0
        else:
            # Otherwise, normalize by normalizing factor
            norm_prob = exp(log_probs[curr_entry] - norm_factor)
        # Record unlogged probability
        norm_probs[curr_entry] = norm_prob
    return norm_probs
        

cpdef int[:] \
    sample_from_multinomial(double[:] probs,
                            int N,
                            int[:] samples):
    """
    Sample one element from multinomial probabilities vector.

    Assumes that the probabilities sum to 1.

    Parameters:
    -----------

    probs : array, vector of probabilities
    N : int, number of samples to draw
    """
    cdef int num_elts = probs.shape[0]
    # Current random samples
    cdef int random_sample = 0
    # Counters over number of samples and number of
    # elements in probability vector
    cdef int curr_sample = 0
    cdef int curr_elt = 0
    cdef double rand_val# = rand() / MY_MAX_INT
    # Get cumulative sum of probability vector
    cdef double[:] cumsum = array_utils.get_double_array(num_elts)
    cumsum = stat_helpers.my_cumsum(probs, cumsum)
    for curr_sample in xrange(N):
        # Draw random number
        rand_val = (rand() % MY_MAX_INT) / MY_MAX_INT
        for curr_elt in xrange(num_elts):
            # If the current cumulative sum is greater than the
            # random number, assign it the index
            if cumsum[curr_elt] >= rand_val:
                random_sample = curr_elt
                break
        samples[curr_sample] = random_sample
    return samples



cpdef int[:] \
  init_assignments(int[:, :] reads,
                   int num_reads,
                   int num_isoforms):
    """
    Initialize assignments of reads to isoforms.

    NOTE: This assumes that the reads that are NOT consistent
    with all the isoforms have been thrown out. If they are
    present, they will be skipped
    """
    # Assignments array to return
    cdef int[:] assignments = array_utils.get_int_array(num_reads)
    cdef int curr_read = 0
    cdef int curr_iso = 0
    for curr_read in xrange(num_reads):
        for curr_iso in xrange(num_isoforms):
            # Assign read to the first isoform it is consistent
            # with. 
            if reads[curr_read, curr_iso] == 1:
                assignments[curr_read] = curr_iso
                break
    return assignments


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





