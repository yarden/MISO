import os
import sys
import time

import numpy as np

import scipy
import scipy.misc 
from scipy.special import gammaln

import misopy
import misopy.miso_scores_single as scores_single
import misopy.miso_scores_paired as scores_paired
import misopy.stat_helpers as stat_helpers

# Global data
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
total_reads = num_reads
num_calls = 2000


def dirichlet_lnpdf(alpha, x):
    """
    Substitute for dirichlet_lnpdf of pygsl.
    """
    dir_log_pdf = \
        gammaln(np.sum(alpha)) - sum(gammaln(alpha)) + np.dot((alpha - 1).T, np.log(x).T)
    return dir_log_pdf


def get_reads(num_reads):
    return reads[0:num_reads]


def get_iso_nums(num_reads):
    return isoform_nums[0:num_reads]


def profile_lndirichlet():
    psi_vector = np.array([0.5, 0.5])
    test_array = np.array([1,2,3,4], dtype=np.float)
    scaled_lens = iso_lens - read_len + 1
    num_calls = 350
    # Get reads and isoform assignments
    num_reads = 500
    reads = get_reads(num_reads)
    iso_nums = get_iso_nums(num_reads)
    # Score dirichlet
    print "Benchmarking lndirichlet functions..."
    print stat_helpers.dirichlet_lnpdf(np.array([1, 1], dtype=np.float),
                                 np.array([0.5, 0.5]))
    print dirichlet_lnpdf(np.array([1, 1]), np.array([0.5, 0.5]))


def profile_cumsum():
    psi_vector = np.array([0.5, 0.5])
    test_array = np.array([1,2,3,4], dtype=np.float)
    scaled_lens = iso_lens - read_len + 1
    num_calls = 350
    # Get reads and isoform assignments
    num_reads = 500
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
        cy_result = stat_helpers.my_cumsum(test_array)
    t2 = time.time()
    print "cy result -> ", cy_result
    for n in range(len(np_result)):
        assert np_result[n] == cy_result[n], \
          "Cumsum does not work."
    print "CYTHON took %.2f seconds" %(t2 - t1)


def profile_sample_reassignments():
    psi_vector = np.array([0.5, 0.5])
    test_array = np.array([1,2,3,4], dtype=np.float)
    scaled_lens = iso_lens - read_len + 1
    num_calls = 350
    # Get reads and isoform assignments
    #num_reads = 2
    #reads = np.array([[1, 0], [0, 1]]) 
    #iso_nums = np.array([0, 1])
    num_reads = 400
    reads = get_reads(num_reads)
    iso_nums = get_iso_nums(num_reads)
    # Score dirichlet
    print "Benchmarking lndirichlet functions..."
    print stat_helpers.dirichlet_lnpdf(np.array([1, 1], dtype=np.float),
                                 np.array([0.5, 0.5]))
    print dirichlet_lnpdf(np.array([1, 1]), np.array([0.5, 0.5]))
    print scores_single.py_sample_from_multinomial(np.array([0.1, 0.3, 0.6]),
                                         100)
    log_psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    log_psi_frag = log_psi_frag - scipy.misc.logsumexp(log_psi_frag)
    log_num_reads_possible_per_iso = np.log(scaled_lens)
    print "Calling sample reassignments %d times " \
          "with %d reads" %(num_calls, num_reads)
    t1 = time.time()
    for n in xrange(num_calls):
        result = scores_single.py_sample_reassignments(reads,
                                                psi_vector,
                                                log_psi_frag,
                                                log_num_reads_possible_per_iso,
                                                scaled_lens,
                                                iso_lens,
                                                num_parts_per_isoform,
                                                iso_nums,
                                                num_reads,
                                                read_len,
                                                overhang_len)
    t2 = time.time()
    print "Took %.2f seconds" %(t2 - t1)
    sys.exit(0)


def profile_sample_from_multinomial():
    """
    Multinomial sampling. This is the major bottleneck
    of the sampling functions.
    """
    print "-" * 20
    p = np.array([0.2, 0.1, 0.5])
    N = len(p)
    num_calls = 1000
    num_reads = 1000
    print "Sampling from multinomial for %d x %d times" %(num_reads,
                                                          num_calls)
    t1 = time.time()
    for n in range(num_reads):
        for x in range(num_calls):
            scores_single.py_sample_from_multinomial(p, N)
    t2 = time.time()
    print "  - Sampling from multinomial took %.2f secs" %(t2 - t1)
    
    

def profile_log_score_reads():
    print "-" * 20
    t1 = time.time()
    psi_vector = np.array([0.5, 0.5])
    scaled_lens = iso_lens - read_len + 1
    num_calls = 3000
    print "Profiling log_score_reads for %d calls..." %(num_calls)
    log_num_reads_possible_per_iso = np.log(scaled_lens)
    log_psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    log_psi_frag = log_psi_frag - scipy.misc.logsumexp(log_psi_frag)
    for n in xrange(num_calls):
        v = scores_single.log_score_reads(reads,
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

    
def profile_sum_log_score_reads():
    print "-" * 20
    t1 = time.time()
    psi_vector = np.array([0.5, 0.5])
    scaled_lens = iso_lens - read_len + 1
    num_calls = 3000
    print "Profiling SUM log score reads (%d reads)" %(num_reads)
    print "Profiling SUM log_score_reads for %d calls..." %(num_calls)
    log_num_reads_possible_per_iso = np.log(scaled_lens)
    log_psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    log_psi_frag = log_psi_frag - scipy.misc.logsumexp(log_psi_frag)
    for n in xrange(num_calls):
        v = scores_single.sum_log_score_reads(reads,
                                       isoform_nums,
                                       num_parts_per_isoform,
                                       iso_lens,
                                       log_num_reads_possible_per_iso,
                                       num_reads,
                                       read_len,
                                       overhang_len)
    t2 = time.time()
    print "SUM log_score_reads took %.2f seconds per %d calls." %(t2 - t1,
                                                                  num_calls)
    

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
    print "-" * 20
    psi_vector = np.array([0.5, 0.5])
    scaled_lens = iso_lens - read_len + 1
    num_calls = 1000
    print "Profiling log score assignments (%d reads)" %(num_reads)
    print "Profiling log score assignments in PYTHON..."
    t1 = time.time()
    for n in range(num_calls):
        v1 = log_score_assignments(isoform_nums,
                                   psi_vector,
                                   scaled_lens,
                                   num_reads)
    t2 = time.time()
    print "Python took %.2f seconds" %(t2 - t1)
    print "Profiling log score assignments in cython..."
    log_psi_frag = np.log(psi_vector) + np.log(scaled_lens)
    log_psi_frag = log_psi_frag - scipy.misc.logsumexp(log_psi_frag)
    t1 = time.time()
    for n in range(num_calls):
        v2 = scores_single.py_log_score_assignments(isoform_nums,
                                             log_psi_frag,
                                             num_reads)
    t2 = time.time()
    print "Cython took %.2f seconds" %(t2 - t1)
    print "RESULTS"
    print "-" * 4
    print "Python: "
    print v1
    print "Cython: "
    print v2
    


def main():
    #profile_sample_from_multinomial()
    profile_sample_reassignments()
    # read scoring
    #profile_log_score_reads()
    #profile_sum_log_score_reads()

    # assignment scoring
    profile_log_score_assignments()

if __name__ == "__main__":
    main()
