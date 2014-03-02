cimport cython
cdef double PI = 3.141592654
cpdef double my_logsumexp(double[:] log_vector,
                         int vector_len)
cpdef double[:] logit(double[:] p,
                      int p_len)
cpdef double[:] logit_inv(double[:] p,
                          int p_len)
