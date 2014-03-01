##
## Test math functions
##
import os
import sys
import time
import unittest

import numpy as np

import scipy
import scipy.misc 
from scipy.special import gammaln

import misopy
import misopy.pyx
import misopy.pyx.matrix_utils as matrix_utils
import misopy.pyx.sampling_utils as sampling_utils


class TestMath(unittest.TestCase):
    """
    Test mathematics functions.
    """
    def setUp(self):
        pass

    
    def test_mat_trans(self):
        """
        Test matrix transpose.
        """
        print "Testing matrix transpose"
        A = np.array([[1,2],
                      [3,4],
                      [5,6]], dtype=float)
        A_trans_pyx = matrix_utils.py_mat_trans(A, A.shape[0], A.shape[1])
        A_trans_numpy = A.T
        assert (np.array_equal(A_trans_pyx, A_trans_numpy)), \
          "Matrix transpose failed (1)."
          
        B = np.array([[1,2,3,4,5],
                      [10,20,30,40,50]], dtype=float)
        B_trans_pyx = matrix_utils.py_mat_trans(B, B.shape[0], B.shape[1])
        B_trans_numpy = B.T
        assert (np.array_equal(B_trans_pyx, B_trans_numpy)), \
          "Matrix transpose failed (2)."


    def test_mat_plus_mat(self):
        print "Testing matrix addition"
        A = np.array([[1, 2, 3],
                      [4, 5, 6]], dtype=float)
        B = np.array([[10, 100, 1000],
                      [1, 1, 1]], dtype=float)
        numpy_C = np.array(np.matrix(A) + np.matrix(B))
        pyx_C = matrix_utils.py_mat_plus_mat(A, A.shape[0], A.shape[1],
                                             B, B.shape[0], B.shape[1])
        assert (np.array_equal(numpy_C, pyx_C)), \
          "Matrix addition failed."

          
    def test_mat_times_mat(self):
        print "Testing matrix multiplication"
        A = np.array([[1, 2, 3],
                      [4, 5, 6]], dtype=float)
        B = np.array([[10, 20],
                      [40, 50],
                      [0, 1]], dtype=float)
        numpy_A_times_B = np.matrix(A) * np.matrix(B)
        pyx_A_times_B = \
          matrix_utils.py_mat_times_mat(A, 
                                        A.shape[0],
                                        A.shape[1],
                                        B.shape[1],
                                        B)
        assert (np.array_equal(pyx_A_times_B, numpy_A_times_B)), \
          "Matrix multiplication failed."
        # Multiply a matrix by a column vector
        A = np.matrix(np.array([[1,2,3],
                                [4,5,6],
                                [7,8,9]]), dtype=float)
        c = np.matrix(np.array([1,2,3], dtype=float)).T
        pyx_A_times_c = A * c
        numpy_A_times_c = \
          matrix_utils.py_mat_times_mat(A, A.shape[0], A.shape[1],
                                        c.shape[1], c)
        assert (np.array_equal(pyx_A_times_c, numpy_A_times_c)), \
          "Matrix multiplication by vector column failed."


    def test_sample_multivar_normal(self):
        mu = np.array([1.1, 1.0], dtype=float)
        sigma = np.matrix(np.array([[0.05, 0],
                                    [0, 0.05]], dtype=float))
        # Get Cholesky decomposition L of Sigma covar matrix
        L = np.linalg.cholesky(sigma)
        print "Cholesky L: "
        print L
        k = mu.shape[0]
        npy_samples = np.random.multivariate_normal(mu, sigma)
        # Cython interface expects mu as a *column* vector
        mu_col = np.matrix(mu).T
        print "mu_col: ", mu_col
        print "Passing L: ", L
        print "mu_col: ", mu_col.shape
        print mu_col.ndim, " <<"
        pyx_samples = sampling_utils.py_sample_multivar_normal(mu_col, L, k)
        print "Numpy samples:"
        print npy_samples
        print "Cython samples:"
        print pyx_samples
        


def main():
    unittest.main()


if __name__ == "__main__":
    main()
