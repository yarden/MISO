cimport numpy as np
ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_float_t
cdef np.float_t sum_array(np.ndarray[np.float_t, ndim=1] input_array, np.int_t array_len)
cdef np.ndarray[double, ndim=1] \
  vect_prod(np.ndarray[double, ndim=1] my_vect,
            int vect_len)
cdef np.ndarray[double, ndim=1] \
  log_vect(np.ndarray[double, ndim=1] my_vect,
           int vect_len)
cdef DTYPE_t array_len(np.ndarray[double, ndim=1] my_array)
cdef np.ndarray[double, ndim=2] \
  mat_times_mat(np.ndarray[double, ndim=2] A,
                int m,
                int n,
                int p,
                np.ndarray[double, ndim=2] B)
cdef np.ndarray[double, ndim=2] \
  mat_dotprod(np.ndarray[double, ndim=2] A,
              int m,
              int n,
              int p,
              np.ndarray[double, ndim=2] B)
cdef np.ndarray[double, ndim=2] \
  mat_plus_mat(np.ndarray[double, ndim=2] A,
               int m,
               int n,
               np.ndarray[double, ndim=2] B,
               int p,
               int q)
cdef np.ndarray[double, ndim=2] \
  row_to_col_vect(np.ndarray[double, ndim=1] row_vect,
                  int k)
cdef np.ndarray[double, ndim=2] mat_trans(np.ndarray[double, ndim=2] A, int m, int n)












