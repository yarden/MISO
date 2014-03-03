##
## Test MISO Engine
##
import os
import sys
import time
import unittest

import numpy as np

import misopy
import misopy.internal_tests
import misopy.internal_tests.test_cases as test_cases

import misopy.pyx
import misopy.pyx.miso_engine as miso_engine


class TestEngine(unittest.TestCase):
    """
    Test MISO scoring functions.
    """
    def setUp(self):
        pass

    def test_engine_single_end(self):
        reads = test_cases.READS
        iso_lens = test_cases.iso_lens
        read_len = 50
        sigma = 0.05
        covar_mat = np.matrix([[sigma]], dtype=float)
        covar_L = np.linalg.cholesky(covar_mat)
        se_engine = miso_engine.SingleEndEngine(reads,
                                                iso_lens,
                                                read_len,
                                                covar_mat,
                                                covar_L,
                                                num_iters=800)
        print "SE engine: "
        print se_engine
        print se_engine.reads
        # Run sampler
        print "Running sampler"
        se_engine.run_sampler()

    def test_engine_paired_end(self):
        pass


def main():
    unittest.main()
    

if __name__ == "__main__":
    main()

  

    

   

    

