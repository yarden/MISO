##
## Define lapack functions to be used
##
cdef extern from "cblas.h":
  enum CBLAS_ORDER:
    CblasRowMajor, CblasColMajor

cdef extern from "f2c.h":
   ctypedef int integer
   ctypedef double doublereal


# Import lapack functions. Note that 'doublereal' type is just 'double'
# when importing.
cdef extern from "clapack.h":
   cdef extern from "f2c.h":
       pass
#   ctypedef int integer
#   ctypedef double doublereal
   integer the_dgemm "dgemm_"(char *transa, char *transb, integer *m, integer *
                              n, integer *k, doublereal *alpha, doublereal *a, integer *lda,
                              doublereal *b, integer *ldb, doublereal *beta, doublereal *c__,
                              integer *ldc)
