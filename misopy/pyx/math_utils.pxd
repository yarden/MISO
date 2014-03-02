cimport cython
cdef double PI = 3.141592654
cpdef double max_val(double m, double n)
cpdef double min_val(double m, double n)
cpdef double my_logsumexp(double[:] log_vector,
                         int vector_len)
cpdef double[:] logit(double[:] p,
                      int p_len)
cpdef double[:] logit_inv(double[:] p,
                          int p_len)
