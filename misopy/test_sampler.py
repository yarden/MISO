import unittest

import scipy
import misopy
import numpy
import numpy as np
import misopy.reads_utils as reads_utils
from misopy.read_simulator import simulate_reads, print_reads_summary, \
                                  read_counts_to_read_list, \
                                  get_reads_summary
from scipy.misc import logsumexp
from scipy.special import gammaln
import misopy.hypothesis_test as ht
from misopy.Gene import Gene, Exon
from misopy.py2c_gene import *

# C MISO interface
import pysplicing

# Cython interface
import misopy.miso_scores_single as scores_single
import misopy.miso_scores_paired as scores_paired

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


class TestSampler(unittest.TestCase):
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
        self.log_psi_frag = \
          np.log(self.psi_vector) + np.log(self.scaled_lens)
        self.log_psi_frag = \
          self.log_psi_frag - scipy.misc.logsumexp(self.log_psi_frag)
        
        

