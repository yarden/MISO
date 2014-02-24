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
