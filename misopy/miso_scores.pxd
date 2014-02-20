cimport libc.math

cdef extern from "math.h":
    double NAN
    double INFINITY


cdef double dirichlet_log_pdf_raw(
    int D,
    double* alpha, int alpha_stride,
    double* vector, int vector_stride,
    )
