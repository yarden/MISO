##
## Test scoring functions
##
import os
import sys
import time
import unittest

import numpy as np

import scipy
import scipy.misc 
from scipy.special import gammaln

import miso_scores_single as scores_single
import miso_scores_paired as scores_paired

num_inc = 3245
num_exc = 22
num_com = 39874
READS = [[1,0]] * num_inc + \
        [[0,1]] * num_exc + \
        [[1,1]] * num_com
READS = np.array(READS, dtype=np.int)
read_len = 40
overhang_len = 4
num_parts_per_iso = np.array([3, 2], dtype=np.int)
iso_lens = np.array([1253, 1172], dtype=np.int)
# Assignment of reads to isoforms: assign half of
# the common reads to isoform 0, half to isoform 1
iso_nums = [0]*3245 + [1]*22 + [0]*19937 + [1]*19937
iso_nums = np.array(iso_nums, dtype=np.int)
num_reads = len(READS)


class TestScores(unittest.TestCase):
    """
    Test MISO scoring functions.
    """
    def setUp(self):
        self.reads = READS
        self.read_len = read_len
        self.overhang_len = overhang_len
        self.num_parts_per_iso = num_parts_per_iso
        self.iso_lens = iso_lens
        self.scaled_lens = self.iso_lens - read_len + 1
        self.log_num_reads_possible_per_iso = np.log(self.scaled_lens)
        self.iso_nums = iso_nums
        self.num_reads = len(self.reads)
        self.psi_vector = np.array([0.8, 0.2])
        # Compute log psi frag
        self.log_psi_frag = np.log(self.psi_vector) + np.log(self.scaled_lens)
        self.log_psi_frag = self.log_psi_frag - scipy.misc.logsumexp(self.log_psi_frag)


    def ____test_log_score_reads(self):
        # Take the first two reads
        curr_num_reads = 2
        two_reads = self.reads[0:2]
        # Check identity of reads
        assert(np.array_equal(two_reads[0], np.array([1, 0])))
        assert(np.array_equal(two_reads[1], np.array([1, 0])))
        # Score the reads given an isoform assignment
        total_log_read_prob = \
          scores_single.sum_log_score_reads(two_reads,
                                     iso_nums[0:2],
                                     num_parts_per_iso,
                                     self.iso_lens,
                                     self.log_num_reads_possible_per_iso,
                                     curr_num_reads,
                                     self.read_len,
                                     self.overhang_len)
        print "TOTAL LOG READ PROB: "
        print total_log_read_prob
        # Compute it by hand: probability of a read is 1 / (# possible positions)
        log_prob_read_1 = np.log(1 / float(self.scaled_lens[iso_nums[0]]))
        log_prob_read_2 = np.log(1 / float(self.scaled_lens[iso_nums[1]]))
        print "log_prob_read_1: ", log_prob_read_1
        print "log_prob_read_2: ", log_prob_read_2
        print log_prob_read_1 + log_prob_read_2
        assert (total_log_read_prob == (log_prob_read_1 + log_prob_read_2)), \
           "Failed to score reads correctly."


    def approx_eq(self, p1, p2, error=0.0001):
        return (np.abs(p1 - p2) < error)


    def __test_log_score_assignments(self):
        curr_num_reads = 2
        two_reads = self.reads[0:curr_num_reads]
        assert self.approx_eq(sum(psi_frag), 1.0), "Psi frag does not sum to 1."
        print "ISO NUMS: ", self.iso_nums[0:curr_num_reads]
        print "LOG PSI FRAG: ", self.log_psi_frag
        print np.log(psi_frag)
        assert (self.approx_eq(self.log_psi_frag[0], np.log(psi_frag)[0])), \
          "Log psi frag not set properly."
        total_log_assignments_prob = \
          scores_single.sum_log_score_assignments(self.iso_nums[0:curr_num_reads],
                                           self.log_psi_frag,
                                           curr_num_reads)
        print "TOTAL LOG ASSIGNMENTS PROB: "
        print total_log_assignments_prob
        # Compute the probability of assignments
        manual_result = (np.log(psi_frag[self.iso_nums[0]]) + \
                         np.log(psi_frag[self.iso_nums[1]]))
        assert (self.approx_eq(manual_result, total_log_assignments_prob)), \
          "Failed to score assignments correctly."


    def test_my_logsumexp(self):
        vals_to_test = [np.array([-1462.26, -1 * np.inf]),
                        np.array([0.1, 0.5])]
        for v in vals_to_test:
            scipy_result = scipy.misc.logsumexp(v)
            result = scores_single.py_my_logsumexp(v, len(v))
            assert (self.approx_eq(scipy_result, result)), \
              "My logsumexp failed on %s" %(str(v))


    def test_sample_reassignment(self):
        curr_num_reads = 200
        subset_reads = self.reads[0:curr_num_reads]
        psi_frag_numer = \
          np.array([(self.scaled_lens[0] * self.psi_vector[0]),
                    (self.scaled_lens[1] * self.psi_vector[1])])
        psi_frag_denom = np.sum(psi_frag_numer)
        psi_frag = psi_frag_numer / psi_frag_denom
        log_psi_frag = np.log(psi_frag)
        result = scores_single.py_sample_reassignments(subset_reads,
                                                       self.psi_vector,
                                                       log_psi_frag,
                                                       self.log_num_reads_possible_per_iso,
                                                       self.scaled_lens,
                                                       self.iso_lens,
                                                       self.num_parts_per_iso,
                                                       self.iso_nums[0:curr_num_reads],
                                                       curr_num_reads,
                                                       self.read_len,
                                                       self.overhang_len)


    def test_init_assignments(self):
        reads = self.reads
        assignments = scores_single.py_init_assignments()
        print "ASSIGNMENTS: ", assignments
        
        


def main():
    unittest.main()


if __name__ == "__main__":
    main()

