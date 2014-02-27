##
## Utilities for working with matrices
##
## Yarden Katz <yarden@mit.edu>
##
cimport numpy as np
import numpy as np

##
## Matrix multiplication from C++
##
#cdef extern from "matrix.h" namespace "matrix":
#   void matrix_mult(double *A, int m, int n, int p, double *B, double *C)
#   void test_mat(double *A, int m, int n)

##
## Matrix multiplication
##
#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef np.ndarray[double, ndim=2] \
  mat_times_mat(np.ndarray[double, ndim=2] A,
                int m,
                int n,
                int p,
                np.ndarray[double, ndim=2] B):
    """
    Matrix x matrix multiplication.

    A : (n x m) matrix
    B : (m x p) matrix

    return C, an (n x p) matrix.
    """
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    # Result matrix
    cdef np.ndarray[double, ndim=2] C = \
      np.zeros((n, p), dtype=float)
    for i in xrange(m):
        for j in xrange(p):
            for k in xrange(n):
                C[i, j] += A[i, k] * B[k, j]
    return C

def py_mat_times_mat(np.ndarray[double, ndim=2] A,
                     int m,
                     int n,
                     int p,
                     np.ndarray[double, ndim=2] B):
    """
    Python interface to mat_times_mat.
    """
    return mat_times_mat(A, m, n, p, B)

##
## Matrix transpose
##
cdef np.ndarray[double, ndim=2] \
  mat_trans(np.ndarray[double, ndim=2] A,
            int m,
            int n):
    """
    Matrix transpose.
 
    A : (m x n) matrix

    Returns A'.
    """
    cdef int i = 0
    cdef int j = 0
    cdef np.ndarray[double, ndim=2] A_trans = np.empty((n, m), dtype=float)
    # By row
    for i in xrange(m):
        # By column
        for j in xrange(n):
            A_trans[j, i] = A[i, j]
    return A_trans

def py_mat_trans(np.ndarray[double, ndim=2] A,
                 int m,
                 int n):
    """
    Python interface to mat_trans.
    """
    return mat_trans(A, m, n)
      
