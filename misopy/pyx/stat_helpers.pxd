cimport libc.math

cdef extern from "math.h":
    double NAN
    double INFINITY


cdef double dirichlet_log_pdf_raw(
    int D,
    double* alpha, int alpha_stride,
    double* vector, int vector_stride,
    )

cpdef double dirichlet_lnpdf(double[:] alpha, double[:] vector)
cpdef double[:] my_cumsum(double[:] input_array, double[:] cumsum_array)
