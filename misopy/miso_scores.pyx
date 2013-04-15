##
## MISO scoring functions in Cython for MCMC sampler
##
import numpy as np
cimport numpy as np

def log_score_reads(np.array reads, np.array isoform_nums,
                    np.array num_parts_per_isoform, np.array iso_lens):
    """
    Score a set of reads given their isoform assignments.
    """
    # The probability of a read being assigned to an isoform that
    # could not have generated it (i.e. where the read is not a
    # substring of the isoform) is zero.  Check for consistency
    overhang_excluded = \
        2*(self.overhang_len - 1)*(num_parts_per_isoform[isoform_nums] - 1)
    # The number of reads possible is the number of ways a 36 nt long read can be
    # aligned onto the length of the current isoforms, minus the number of positions
    # that are ruled out due to overhang constraints.
    num_reads_possible = \
        (iso_lens[isoform_nums] - self.read_len + 1) - overhang_excluded
    log_prob_reads = log(1) - log(num_reads_possible)
    zero_prob_indx = nonzero(reads[xrange(self.num_reads), isoform_nums] == 0)[0]
    # Assign probability 0 to reads inconsistent with assignment
    log_prob_reads[zero_prob_indx] = -inf
    return log_prob_reads


def primes(int kmax):
    cdef int n, k, i
    cdef int p[1000]
    result = []
    if kmax > 1000:
        kmax = 1000
    k = 0
    n = 2
    while k < kmax:
        i = 0
        while i < k and n % p[i] != 0:
            i = i + 1
        if i == k:
            p[k] = n
            k = k + 1
            result.append(n)
        n = n + 1
    return result
