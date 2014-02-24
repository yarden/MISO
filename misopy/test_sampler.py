import unittest

import os
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
import misopy.Gene as gene_utils

import misopy.miso_sampler as miso_sampler
import misopy.settings as settings

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
iso_lens = np.array([1200, 1000], dtype=np.int)
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
        self.num_isoforms = len(self.iso_lens)
        self.psi_vector = np.array([0.8, 0.2])
        # Compute log psi frag
        self.log_psi_frag = \
          np.log(self.psi_vector) + np.log(self.scaled_lens)
        self.log_psi_frag = \
          self.log_psi_frag - scipy.misc.logsumexp(self.log_psi_frag)
        # Read the params from a file
        self.params = None
        self.settings_fname = \
          os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "settings",
                       "miso_settings.txt")
        print "Settings filename: %s" %(self.settings_fname)
        self.settings = settings.load_settings(self.settings_fname)
        self.params = settings.Settings.get_sampler_params()
        print self.settings
        print self.params
        self.exon_lens = [500, 200, 500]
        self.isoforms = [[1, 2, 3],
                         [1, 3]]
#        self.gene = gene_utils.load_genes_from_gff()


    def test_sampler1(self):
        sampler_obj = miso_sampler.MISOSampler(self.params,
                                               paired_end=False,
                                               log_dir="./miso_logs/")
        print "sampler_obj: ", sampler_obj
        output_fname = "./sampler_output"
        num_iters = 1000
        reads = self.reads
        gene = self.gene
        hyperparameters = np.ones(self.num_isoforms)
        sampler_obj.run_sampler(self,
                                num_iters,
                                reads,
                                gene,
                                hyperparameters,
                                self.params,
                                output_fname)



def main():
    unittest.main()
        

if __name__ == "__main__":
    main()
        

