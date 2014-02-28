##
## Utilities for working with matrices
##
## Yarden Katz <yarden@mit.edu>
##
import numpy as np
cimport numpy as np

from libc.math cimport log
from libc.math cimport exp

cimport lapack

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_float_t


##
## Vector utilities
##
cdef np.ndarray[double, ndim=1] \
  vect_prod(np.ndarray[double, ndim=1] my_vect,
            int vect_len):
    """
    Return vector product.
    """
    cdef int i = 0
    cdef double prod_result = 1.0
    for i in xrange(my_vect):
        prod_result = prod_result * my_vect[i]
    return prod_result

def my_vect_prod(np.ndarray[double, ndim=1] my_vect,
                 int vect_len):
    return vect_prod(my_vect, vect_len)


cdef DTYPE_float_t \
  sum_array(np.ndarray[DTYPE_float_t, ndim=1] input_array,
            DTYPE_t array_len):
    cdef DTYPE_t j = 0
    cdef DTYPE_float_t result = 0.0
    for j in xrange(array_len):
        result += input_array[j]
    return result


cdef np.ndarray[double, ndim=1] \
  log_vect(np.ndarray[double, ndim=1] my_vect,
           int vect_len):
    """
    Return log of vector
    """
    cdef int i = 0
    cdef np.ndarray[double, ndim=1] log_my_vect = \
      np.empty(vect_len, dtype=float)
    for i in xrange(vect_len):
        log_my_vect[i] = log(my_vect[i])
    return log_my_vect


cdef DTYPE_t array_len(np.ndarray[double, ndim=1] my_array):
    """
    Return length of 1d array.
    """
    return my_array.shape[0]

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
## Matrix dot product (for 2d arrays, this is equivalent to
## matrix multiplication.)
##
cdef np.ndarray[double, ndim=2] \
  mat_dotprod(np.ndarray[double, ndim=2] A,
              int m,
              int n,
              int p,
              np.ndarray[double, ndim=2] B):
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

##
## Cholesky decomposition 
##
# cdef np.ndarray[double, ndim=2] \
#   la_cholesky(np.ndarray[double, ndim=2] A):
#     """
#     Cholesky decomposition using CLAPACK.

#     Given input matrix A (m x n), find:
    
#     A = LL'

#     Return L. Assumes A and L are 2d double arrays.

#     NOTE: Assumes A is a C array (not numpy array!)
#     but returns a numpy array in result.
#     """
#     # CLAPACK takes "integer" types (i.e. a "long int")
#     cdef long int m = <long int>A.shape[0]
#     cdef long int n = <long int>A.shape[1]
#     cdef np.ndarray[double, ndim=2, mode="c"] L = \
#       np.empty((A.shape[0], A.shape[1]), dtype=float)
#     return ;
      
