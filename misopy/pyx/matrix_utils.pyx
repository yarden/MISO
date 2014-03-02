##
## Utilities for working with matrices
##
## Yarden Katz <yarden@mit.edu>
##
from libc.math cimport log
from libc.math cimport exp

cimport lapack

##
## Vector utilities
##
cpdef double \
  vect_prod(double[:] my_vect,
            int vect_len):
    """
    Return vector product.
    """
    cdef int i = 0
    cdef double prod_result = 1.0
    for i in xrange(my_vect.shape[0]):
        prod_result = prod_result * my_vect[i]
    return prod_result


cpdef double \
  sum_array(double[:] input_array,
            int array_len):
    cdef int j = 0
    cdef double result = 0.0
    for j in xrange(array_len):
        result += input_array[j]
    return result


cpdef double[:] \
  log_vect(double[:] my_vect,
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


cpdef int array_len(double[:] my_array):
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
## Matrix addition
##
cpdef double[:, :] \
  mat_plus_mat(double[:, :] A,
               int m,
               int n,
               double[:, :] B,
               int p,
               int q):
    """
    Add two matrices together. Adds matrix A (m x n)
    with matrix B (p x q).
    """
    cdef np.ndarray[double, ndim=2] added_mat = \
      np.empty((m, n), dtype=float)
    cdef int i = 0
    cdef int j = 0
    for i in xrange(m):
        for j in xrange(n):
            added_mat[i][j] = A[i][j] + B[i][j]
    return added_mat


cpdef double[:, :] \
  row_to_col_vect(double[:] row_vect,
                  int k):
    """
    Convert row vector (1d) to column vector (2d).

    This takes a 1d np.ndarray of length k and converts it
    to a 2d matrix of (k x 1).
    """
    cdef int i = 0
    cdef np.ndarray[double, ndim=2] col_vect = \
      np.empty((k, 1), dtype=float)
    for i in xrange(k):
        col_vect[i][0] = row_vect[i]
    return col_vect
    

##
## Matrix multiplication
##
#@cython.boundscheck(False)
#@cython.wraparound(False)
cpdef double[:, :] \
  mat_times_mat(double[:, :] A,
                int m,
                int n,
                int p,
                double[:, :] B):
    """
    Matrix x matrix multiplication.

    A : (m x n) matriax
    B : (n x p) matrix

    return C, an (n x p) matrix.
    """
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    # Result matrix
    cdef np.ndarray[double, ndim=2] C = \
      np.zeros((m, p), dtype=float)
    for i in xrange(m):
        for j in xrange(n):
            for k in xrange(p):
                C[i, k] += A[i, j] * B[j, k]
    return C


##
## Matrix dot product (for 2d arrays, this is equivalent to
## matrix multiplication.)
##
cpdef double[:, :] \
  mat_dotprod(double[:, :] A,
              int m,
              int n,
              int p,
              double[:, :] B):
    """
    Dot product of matrix A x matrix B.
    """
    return mat_times_mat(A, m, n, p, B)

##
## Matrix times column vector
##
# ...


##
## Matrix transpose
##
cpdef double[:, :] \
  mat_trans(double[:, :] A,
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

      
