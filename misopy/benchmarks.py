import os
import sys
import time

import numpy as np

import scipy
import scipy.misc 
from scipy.special import gammaln

import miso_scores as scores

num_inc = 3245
num_exc = 22
num_com = 39874
reads = [[1,0]] * num_inc + \
        [[0,1]] * num_exc + \
        [[1,1]] * num_com
reads = np.array(reads, dtype=np.int)
isoform_nums = []
read_len = 40
overhang_len = 4
num_parts_per_isoform = np.array([3, 2], dtype=np.int)
iso_lens = np.array([1253, 1172], dtype=np.int)
# Assignment of reads to isoforms: assign half of
# the common reads to isoform 0, half to isoform 1
isoform_nums = [0]*3245 + [1]*22 + [0]*19937 + [1]*19937
isoform_nums = np.array(isoform_nums, dtype=np.int)
num_reads = len(reads)
num_calls = 2000


def dirichlet_lnpdf(alpha, x):
    """
    Substitute for dirichlet_lnpdf of pygsl.
    """
    dir_log_pdf = \
        gammaln(np.sum(alpha)) - sum(gammaln(alpha)) + np.dot((alpha - 1).T, np.log(x).T)
    return dir_log_pdf


def profile_sample_reassignments():
    psi_vector = np.array([0.5, 0.5])
    test_array = np.array([1,2,3,4], dtype=np.float)
    scaled_lens = iso_lens - read_len + 1
    iso_nums = isoform_nums
    print "Profiling numpy cumsum"
    t1 = time.time()
    np_result = None
    for n in np.arange(num_calls):
        np_result = np.cumsum(test_array)
    t2 = time.time()
    print "Took %.2f seconds" %(t2 - t1)
    print "np result -> ", np_result
    print "Profiling CYTHON cumsum"
    t1 = time.time()
    cy_result = None
    for n in np.arange(num_calls):
        cy_result = scores.my_cumsum(test_array)
    t2 = time.time()
    print "cy result -> ", cy_result
    for n in range(len(np_result)):
        assert np_result[n] == cy_result[n], \
          "Cumsum does not work."
    print "CYTHON took %.2f seconds" %(t2 - t1)
    # Score dirichlet
    print "Benchmarking lndirichlet functions..."
    print scores.dirichlet_lnpdf(np.array([1, 1], dtype=np.float),
                                 np.array([0.5, 0.5]))
    print dirichlet_lnpdf(np.array([1, 1]), np.array([0.5, 0.5]))
    print scores.py_sample_from_multinomial(np.array([0.1, 0.3, 0.6]),
                                         100)
    print "psi vector: ", psi_vector
    print "iso_lens: ", iso_lens
    print "scaled_lens: ", scaled_lens
    print "num parts: ", num_parts_per_isoform
    log_psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    log_psi_frag = log_psi_frag - scipy.misc.logsumexp(log_psi_frag)
    log_num_reads_possible_per_iso = np.log(scaled_lens)
    t1 = time.time()
    result = scores.sample_reassignments(reads,
                                         psi_vector,
                                         log_psi_frag,
                                         log_num_reads_possible_per_iso,
                                         scaled_lens,
                                         iso_lens,
                                         num_parts_per_isoform,
                                         isoform_nums,
                                         num_reads,
                                         read_len,
                                         overhang_len)
    t2 = time.time()
    print "Took %.2f seconds" %(t2 - t1)
    sys.exit(0)


def profile_sample_from_multinomial():
    p = np.array([0.2, 0.1, 0.5])
    N = len(p)
    print "Sampling from multinomial for %d times" %(num_reads)
    t1 = time.time()
    for n in range(num_reads):
        scores.py_sample_from_multinomial(p, N)
    t2 = time.time()
    print "Sampling from multinomial took %.2f secs" %(t2 - t1)
    
    

def profile_log_score_reads():
    t1 = time.time()
    psi_vector = np.array([0.5, 0.5])
    scaled_lens = iso_lens - read_len + 1
    num_calls = num_reads
    print "Profiling log_score_reads for %d calls..." %(num_calls)
    log_num_reads_possible_per_iso = np.log(scaled_lens)
    for n in xrange(num_calls):
        ###
        ### TODO: TEST LOG SCORE READS HERE TO MAKE SURE ISO NUMS WORKS
        ###
        log_psi_frag = np.log(psi_vector) + np.log(scaled_lens)
        log_psi_frag = log_psi_frag - scipy.misc.logsumexp(log_psi_frag)
        v = scores.log_score_reads(reads,
                                  isoform_nums,
                                  num_parts_per_isoform,
                                  iso_lens,
                                  log_num_reads_possible_per_iso,
                                  num_reads,
                                  read_len,
                                  overhang_len)
    t2 = time.time()
    print "log_score_reads took %.2f seconds per %d calls." %(t2 - t1,
                                                              num_calls)
    

def py_log_score_reads(reads,
                       isoform_nums,
                       num_parts_per_isoform,
                       iso_lens,
                       read_len,
                       overhang_len,
                       num_reads):
    """
    Score a set of reads given their isoform assignments.
    """
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
    

def log_score_assignments(isoform_nums, psi_vector, scaled_lens, num_reads):
    """
    Score an assignment of a set of reads given psi
    and a gene (i.e. a set of isoforms).
    """
    psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    psi_frag = psi_frag - scipy.misc.logsumexp(psi_frag)
    psi_frags = np.tile(psi_frag, [num_reads, 1])
    return psi_frags[np.arange(num_reads), isoform_nums]


def profile_log_score_assignments():
    psi_vector = np.array([0.5, 0.5])
    scaled_lens = iso_lens - read_len + 1
    num_calls = 1000
    print "Profiling log score assignments in PYTHON..."
    t1 = time.time()
    for n in range(num_calls):
        v1 = log_score_assignments(isoform_nums, psi_vector, scaled_lens, num_reads)
    t2 = time.time()
    print "Took %.2f seconds" %(t2 - t1)
    print "Profiling log score assignments in cython..."
    t1 = time.time()
    for n in range(num_calls):
        v2 = scores.log_score_assignments(isoform_nums, psi_vector, scaled_lens, num_reads)
    t2 = time.time()
    print "Took %.2f seconds" %(t2 - t1)
    print "Profiling LOOP log score assignments in cython..."
    t1 = time.time()
    # Precompute psi_frag
    log_psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    log_psi_frag = log_psi_frag - scipy.misc.logsumexp(log_psi_frag)
    for n in range(num_calls):
        v3 = scores.loop_log_score_assignments(isoform_nums, log_psi_frag, num_reads)
    t2 = time.time()
    print "LOOP Took %.2f seconds" %(t2 - t1)
    print "RESULTS"
    print "-" * 4
    print v1, v2, v3
    


def main():
    profile_sample_from_multinomial()
    profile_sample_reassignments()
    profile_log_score_reads()
    #profile_log_score_assignments()

if __name__ == "__main__":
    main()
