##
## Profile the MISO scoring functions
##
import os
import sys
import time

import numpy as np

import misopy
import misopy.miso_scores as miso_scores


def profile_log_score_reads():
    t1 = time.time()
    num_inc = 3245
    num_exc = 22
    num_com = 39874
    reads = [[1,0]] * num_inc + \
            [[0,1]] * num_exc + \
            [[1,1]] * num_com
    reads = np.array(reads)
    isoform_nums = []
    read_len = 40
    overhang_len = 4
    num_parts_per_isoform = np.array([3, 2], dtype=int)
    iso_lens = np.array([1253, 1172], dtype=int)
    # Assignment of reads to isoforms: assign half of
    # the common reads to isoform 0, half to isoform 1
    isoform_nums = [0]*3245 + [1]*22 + [0]*19937 + [1]*19937
    isoform_nums = np.array(isoform_nums, dtype=int)
    num_reads = len(reads)
    num_calls = 5000
    print "Profiling log_score_reads for %d calls..." %(num_calls)
    for n in xrange(num_calls):
        miso_scores.log_score_reads(reads,
                                    isoform_nums,
                                    num_parts_per_isoform,
                                    iso_lens,
                                    read_len,
                                    overhang_len,
                                    num_reads)
    t2 = time.time()
    print "log_score_reads took %.2f seconds per %d calls." %(t2 - t1,
                                                              num_calls)


def main():
    profile_log_score_reads()

if __name__ == "__main__":
    main()
