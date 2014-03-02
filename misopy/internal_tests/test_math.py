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
import misopy.pyx.math_utils as math_utils
import misopy.pyx.matrix_utils as matrix_utils
import misopy.pyx.sampling_utils as sampling_utils
import misopy.internal_tests.py_scores as py_scores


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
        A_trans_pyx = np.empty((A.shape[1], A.shape[0]), dtype=float)
        A_trans_pyx = matrix_utils.mat_trans(A, A.shape[0], A.shape[1], A_trans_pyx)
        A_trans_numpy = A.T
        assert (np.array_equal(A_trans_pyx, A_trans_numpy)), \
          "Matrix transpose failed (1)."
          
        B = np.array([[1,2,3,4,5],
                      [10,20,30,40,50]], dtype=float)
        B_trans_pyx = np.empty((B.shape[1], B.shape[0]), dtype=float)
        B_trans_pyx = matrix_utils.mat_trans(B, B.shape[0], B.shape[1],
                                             B_trans_pyx)
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
        pyx_C = np.empty((A.shape[0], A.shape[1]), dtype=float)
        pyx_C = matrix_utils.mat_plus_mat(A, A.shape[0], A.shape[1],
                                          B, B.shape[0], B.shape[1],
                                          pyx_C)
        assert (np.array_equal(numpy_C, pyx_C)), \
          "Matrix addition failed."

          
    def test_mat_times_mat(self):
        print "Testing matrix multiplication"
        A = np.array([[1, 2, 3],
                      [4, 5, 6]], dtype=float)
        B = np.array([[10, 20],
                      [40, 50],
                      [0, 1]], dtype=float)
        C = np.empty((A.shape[0], B.shape[1]), dtype=float)
        numpy_A_times_B = np.matrix(A) * np.matrix(B)
        pyx_A_times_B = \
          matrix_utils.mat_times_mat(A,
                                     A.shape[0],
                                     A.shape[1],
                                     B.shape[1],
                                     B,
                                     C)
        pyx_A_times_B = np.asarray(pyx_A_times_B)
        assert (np.array_equal(pyx_A_times_B, numpy_A_times_B)), \
          "Matrix multiplication failed."
        # Multiply a matrix by a column vector
        A = np.matrix(np.array([[1,2,3],
                                [4,5,6],
                                [7,8,9]]), dtype=float)
        c = np.matrix(np.array([1,2,3], dtype=float)).T
        pyx_A_times_c = np.empty((A.shape[0], 1), dtype=float)
        numpy_A_times_c = A * c
        pyx_A_times_c = \
          matrix_utils.mat_times_mat(A, A.shape[0], A.shape[1],
                                     1, c,
                                     pyx_A_times_c)
        pyx_A_times_c = np.asarray(pyx_A_times_c)
        assert (np.array_equal(pyx_A_times_c, numpy_A_times_c)), \
          "Matrix multiplication by vector column failed."


    def test_sample_multivar_normal(self):
        print "Testing sampling from multivariate normal"
        mu = np.array([2.05, 0.55], dtype=float)
        sigma = np.matrix(np.array([[0.05, 0],
                                    [0, 0.05]], dtype=float))
        # Get Cholesky decomposition L of Sigma covar matrix
        L = np.linalg.cholesky(sigma)
        k = mu.shape[0]
        all_numpy_samples = []
        all_pyx_samples = []
        # Compile a list of all the samples
        num_iter = 1000
        for n in range(num_iter):
            npy_samples = np.random.multivariate_normal(mu, sigma)
            all_numpy_samples.append(npy_samples)
            # Cython interface expects mu as a *column* vector
            mu_col = np.matrix(mu).T
            pyx_samples = sampling_utils.sample_multivar_normal(mu_col, L, k)
            pyx_samples = list(np.asarray(pyx_samples))
            all_pyx_samples.append(pyx_samples)
        # The means should equal the mean we started with
        all_numpy_samples = np.array(all_numpy_samples)
        all_pyx_samples = np.array(all_pyx_samples)
        print "Numpy mean across %d iterations" %(num_iter)
        numpy_mean = np.mean(all_numpy_samples, axis=0)
        print numpy_mean
        print "Cython mean across %d iterations" %(num_iter)
        pyx_mean = np.mean(all_pyx_samples, axis=0)
        print pyx_mean
        # The two methods should yield very similar means
        error_range = 0.025
        assert (py_scores.approx_eq(numpy_mean[0],
                                    pyx_mean[0],
                                    error=error_range) and \
                py_scores.approx_eq(numpy_mean[1],
                                    pyx_mean[1],
                                    error=error_range)), \
                "Numpy and Cython average values of sampled normals are different."


    def test_logit(self):
        print "Testing logit transform"
        x = np.array([0.5, 0.6, 0.01, 0.001, 0.9999, 0, 0.99],
                     dtype=float)
        numpy_logit = py_scores.logit(x)
        pyx_logit = np.asarray(math_utils.logit(x, x.shape[0]))
        print "Numpy logit: "
        print numpy_logit
        print "Pyx logit: "
        print pyx_logit
        assert (np.array_equal(numpy_logit, pyx_logit)), \
          "Logit failed."

          
    def test_logit_inv(self):
        print "Testing inverse logit transform"
        x = np.array([-100, 100, 0.5, 0.6, -0.58, 0.8, 1, 0],
                     dtype=float)
        numpy_logit_inv = py_scores.logit_inv(x)
        pyx_logit_inv = np.asarray(math_utils.logit_inv(x, x.shape[0]))
        print "Numpy logit inv: "
        print numpy_logit_inv
        print "Pyx logit inv: "
        print pyx_logit_inv
        assert (np.array_equal(numpy_logit_inv, pyx_logit_inv)), \
          "Logit inverse failed."
        
        
        

def main():
    unittest.main()


if __name__ == "__main__":
    main()
