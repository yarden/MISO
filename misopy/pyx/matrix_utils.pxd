cimport cython
cpdef double vect_prod(double[:] my_vect, int vect_len)
cpdef double[:] log_vect(double[:] my_vect, int vect_len)
cpdef int array_len(double[:] my_array)
cpdef double[:, :] mat_times_mat(double[:, :] A,
                                 int m,
                                 int n,
                                 int p,
                                 double[:, :] B,
                                 double[:, :] C)
cpdef double[:, :] \
  mat_plus_mat(double[:, :] A,
               int m,
               int n,
               double[:, :] B,
               int p,
               int q,
               double[:, :] added_mat)
cpdef double[:, :] \
  row_to_col_vect(double[:] row_vect,
                  int k)
cpdef double[:, :] mat_trans(double[:, :] A,
                             int m,
                             int n,
                             double[:, :])












