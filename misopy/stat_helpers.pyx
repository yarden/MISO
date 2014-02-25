import numpy as np
cimport numpy as np

np.import_array()
cimport cython

from libc.math cimport log
from libc.math cimport exp
from libc.stdlib cimport rand

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_float_t

cdef float MY_MAX_INT = float(10000)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def my_cumsum(np.ndarray[double, ndim=1] input_array):
    """
    Return cumulative sum of array.
    """
    # Cumulative sum at every position
    cdef np.ndarray[double, ndim=1] cumsum_array = np.empty_like(input_array)
    cdef int curr_elt = 0
    cdef int num_elts = input_array.shape[0]
    # Current cumulative sum: starts at first element
    cdef double curr_cumsum = 0.0
    for curr_elt in xrange(num_elts):
        cumsum_array[curr_elt] = (input_array[curr_elt] + curr_cumsum)
        curr_cumsum = cumsum_array[curr_elt]
    return cumsum_array


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def dirichlet_lnpdf(np.ndarray[double, ndim=1] alpha,
                    np.ndarray[double, ndim=1] vector):
    """
    Wrapper for dirichlet log pdf scoring function.
    """
    cdef int D = vector.size
    return dirichlet_log_pdf_raw(D,
                                 &alpha[0], alpha.strides[0],
                                 &vector[0], vector.strides[0])


# from Borg project
@cython.infer_types(True)
cdef double dirichlet_log_pdf_raw(
    int D,
    double* alpha, int alpha_stride,
    double* vector, int vector_stride,
    ):
    """Compute the log of the Dirichlet PDF evaluated at one vector."""

    cdef void* alpha_p = alpha
    cdef void* vector_p = vector

    # first term
    term_a = 0.0

    for d in xrange(D):
        term_a += (<double*>(alpha_p + alpha_stride * d))[0]

    term_a = libc.math.lgamma(term_a)

    # second term
    term_b = 0.0

    for d in xrange(D):
        term_b += libc.math.lgamma((<double*>(alpha_p + alpha_stride * d))[0])

    # third term
    cdef double alpha_d
    cdef double vector_d

    term_c = 0.0

    for d in xrange(D):
        alpha_d = (<double*>(alpha_p + alpha_stride * d))[0]
        vector_d = (<double*>(vector_p + vector_stride * d))[0]

        term_c += (alpha_d - 1.0) * libc.math.log(vector_d)

    # ...
    return term_a - term_b + term_c


##
## Matrix operations
##
from libc.math cimport pow

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double determinant(double[:, :] M):
    '''
    Compute the determinant of a square nxn matrix using
    the recursive formula:
 
        det(M) = sum_{i=0}^{n-1} (-1^i)*M[0,i]*det(M_minor(0,i))

    where M_minor(0,i) is the (n-1)x(n-1) matrix formed by
    removing row 0 and column i from M.
    
    source: http://nbviewer.ipython.org/github/carljv/cython_testing/blob/master/cython_linalg.ipynb
    '''
    assert M.shape[0] == M.shape[1], 'Matrix is not square.'
    
    cdef int i, j
    cdef int n = M.shape[0]
    cdef double det = 0.0
    cdef double coef
    cdef double[:, :] M_minor = np.empty((n-1, n-1))
    
    if n == 1:
        # If M is a scalar (1x1) just return it
        return M[0, 0]
    else:
        # If M is nxn, then get its (n-1)x(n-1) minors
        # (one for each of M's n columns) and compute 
        # their determinants and add them to the summation.
        for j in xrange(n):
            coef = pow(-1, j) * M[0, j]
            _get_minor(M, M_minor, 0, j)
            det += coef * determinant(M_minor)
        return det

    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _get_minor(double[:, :] M, double[:, :] M_minor, 
                    int row, int col):
    '''
    Return the minor of a matrix, by removing a specified
    row and column

    If M is nxn, then _get_minor(M, row, col) will fill in
    the (n-1)x(n-1) matrix M_minor by removing row `row` 
    and column `col` from M.
    
    source: http://nbviewer.ipython.org/github/carljv/cython_testing/blob/master/cython_linalg.ipynb
    '''
    cdef:
        int n = M.shape[0]
        int i_to, j_to, i_from, j_from
  
        
    # _from indicates the index of the original
    # matrix M, _to, indicates the index of the
    # result matrix M_minor.
    i_from = 0
    for i_to in xrange(n-1):
        if (i_to == row): 
            # This is the row to exclude from the
            # minor, so skip it.
            i_from += 1
        j_from = 0
        for j_to in xrange(n-1):
            if (j_to == col): 
            # This is the column to exclude from the
            # minor, so skip it.
                j_from += 1
            
            M_minor[i_to, j_to] = M[i_from, j_from]
            j_from += 1
        
        i_from += 1
