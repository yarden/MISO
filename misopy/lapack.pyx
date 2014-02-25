##
## Define lapack functions to be used
##
cdef extern from "cblas.h":
  enum CBLAS_ORDER:
    CblasRowMajor, CblasColMajor

    
# Import lapack functions. Note that 'doublereal' type is just 'double'
# when importing.
cdef extern from "clapack.h":
    # int dgemm_(char *transa, char *transb, integer *m, integer *
	# n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
	# doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
	# integer *ldc)
    int C_dgemm "dgemm_"(char *transa,
                         char *transb,
                         int *m,
                         int * n,
                         int *k,
                         double *alpha,
                         double *a,
                         int *lda,
                         double *b,
                         int *ldb,
                         double *beta,
                         double *c__,
                         int *ldc)
print "C_dgemm:"
