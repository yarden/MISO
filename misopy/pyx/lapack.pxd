cimport numpy as np
cdef np.ndarray[double, ndim=2] \
  la_cholesky_decomp(np.ndarray[double, ndim=2] A,
                     int num_rows,
                     int num_cols)
