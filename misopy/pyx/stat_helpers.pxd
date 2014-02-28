cimport numpy as np
cimport libc.math

cdef extern from "math.h":
    double NAN
    double INFINITY


cdef double dirichlet_log_pdf_raw(
    int D,
    double* alpha, int alpha_stride,
    double* vector, int vector_stride,
    )

cdef double dirichlet_lnpdf(np.ndarray[double, ndim=1] alpha,
                            np.ndarray[double, ndim=1] vector)
cdef np.ndarray[double, ndim=1] my_cumsum(np.ndarray[double, ndim=1] input_array)
