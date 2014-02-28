##
## Test calling lapack functions
##
import os
import sys
import time
import unittest

import numpy as np

import misopy
import misopy.pyx
import misopy.pyx.lapack as lapack

class TestLAPACK(unittest.TestCase):
    """
    Test LAPACK.
    """
    def setUp(self):
        pass

    
    def test_cholesky(self):
        """
        Test Cholesky decomposition.
        """
        # Create a C array
        A = \
          np.ascontiguousarray(np.array([[0.05, 0, 0],
                                         [0, 0.05, 0],
                                         [0, 0, 0.05]], dtype=float))
        # Verify we get the same answer as numpy
        A_decomp_np = np.linalg.cholesky(A)
        print "Numpy decomp: "
        print A_decomp_np
        A_decomp_lapack = \
          lapack.py_la_cholesky_decomp(A, A.shape[0], A.shape[1])
        print "LAPACK decomp: "
        print A_decomp_lapack
        assert (np.array_equal(A_decomp_lapack, A_decomp_np)), \
          "LAPACK and Numpy give different decompositions."
        

def main():
    unittest.main()
    

if __name__ == "__main__":
    main()
