cimport numpy as np
import numpy as np
##
## Define lapack functions to be used
##
cdef extern from "f2c.h":
   ctypedef int integer
   ctypedef double doublereal

# Import lapack functions.
cdef extern from "clapack.h":
   integer f2c_dgemm(char *transa, char *transb, integer *m,
                     integer *n, integer *k, doublereal *alpha,
                     doublereal *a, integer *lda, doublereal *b,
                     integer *ldb, doublereal *beta, doublereal *c__,
                     integer *ldc)
   int dpotrf_(char *uplo, integer *n, doublereal *a,
               integer * lda, integer *info)

cdef np.ndarray[double, ndim=2] \
  copy_mat(np.ndarray[double, ndim=2] my_mat):
    """
    Make copy of a 2d array.
    """
    cdef np.ndarray[double, ndim=2] my_new_mat = \
      np.empty((my_mat.shape[0], my_mat.shape[1]), dtype=float)
    cdef int i = 0
    cdef int j = 0
    for i in xrange(my_mat.shape[0]):
        for j in xrange(my_mat.shape[1]):
            my_new_mat[i][j] = my_mat[i][j]
    return my_new_mat

   
##
## Cholesky decomposition
##
cdef np.ndarray[double, ndim=2] \
  la_cholesky_decomp(np.ndarray[double, ndim=2] A,
                     int num_rows,
                     int num_cols):
    """
    Interface to CLAPACK Cholesky decomposition function.

    NOTE: This function destructively modifies the
    input array A!
    """
    # Always decompose matrix to A = LL'
    cdef char uplo = 'L'
    # N is number of columns 
    cdef integer n = <integer>num_cols
    # LDA is number of rows
    cdef integer lda = <integer>num_rows
    # Return value of functio
    cdef integer info
    dpotrf_(&uplo, &n, &A[0,0], &lda, &info)
    if info != 0:
        # Something went wrong
        print "Cholesky decomposition in CLAPACK failed!"
        raise Exception, "CLAPACK failure."
    return A

   
cdef int main():
    cdef char transa_val = 'N'
    cdef char *transa = &transa_val

    cdef char transb_val = 'N'
    cdef char *transb = &transb_val

    #cdef integer m = 3
    #cdef integer n = 3
    #cdef integer k = 3

    cdef double alpha = 1.0
    cdef np.ndarray[double, ndim=2, mode="c"] a = \
     np.asarray(np.array([[1, 2, 3],
                          [4, 5, 6],
                          [7, 8, 9]],
                          dtype=float),
                          order="c")

    cdef np.ndarray[double, ndim=2, mode="c"] b = \
     np.asarray(np.array([[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]],
                          dtype=float), order="c")

    cdef np.ndarray[double, ndim=2, mode="c"] c = \
     np.zeros([3, 3], dtype=float)

    # cdef double a[3][3]
    # a[0][0] = 1
    # a[0][1] = 2
    # a[0][2] = 3
    # a[1][0] = 4
    # a[1][1] = 5
    # a[1][2] = 6
    # a[2][0] = 7
    # a[2][1] = 8
    # a[2][2] = 9
    # cdef double b[3][3]
    # b[0][0] = 1
    # b[0][1] = 0
    # b[0][2] = 0
    # b[1][0] = 0
    # b[1][1] = 0
    # b[1][2] = 0
    # b[2][0] = 1
    # b[2][1] = 1
    # b[2][2] = 1
    # cdef doublereal c[3][3]

    #cdef integer lda = 3
    #cdef integer ldb = 3

    #cdef doublereal beta = 0.0

    #cdef integer ldc = 3

    #print "a[0][0]", a[0][0]

    # Pass C arrays
    #f2c_dgemm(&transa_val, &transb_val, &m, &n, &k, &alpha, &a[0,0], &lda, &b[0,0], &ldb, &beta, &c[0,0], &ldc)    #for i in xrange(3):
    #    for j in xrange(3):
    #        print "C(%d,%d) is %.2f" %(i, j, c[i][j])

    # Try cholesky decomposition
    #f2c_dpotrf(char *uplo, integer *n, doublereal *a, integer *
	#lda, integer *info)
    cdef char uplo = 'L'
    # number of rows
    cdef integer lda = 3
    # number of columns (??)
    cdef integer n = 3
    cdef integer info = 0
    print "Calling LAPACK Cholesky..."
    dpotrf_(&uplo, &n, &a[0,0], &lda, &info)
    print a

    return 0

main()

