
#include "splicing.h"
#include "splicing_error.h"

int splicingdgemv_(char *trans, int *m, int *n, double *alpha,
    double *a, int *lda, double *x, int *incx,
    double *beta, double *y, int *incy);

int splicing_dgemv(int transpose, double alpha,
		   const splicing_matrix_t* a, const splicing_vector_t* x,
		   double beta, splicing_vector_t* y) {
  char trans = transpose ? 'T' : 'N';
  int m, n;
  int inc = 1;

  m = splicing_matrix_nrow(a);
  n = splicing_matrix_ncol(a);

  if (splicing_vector_size(x) != (transpose ? m : n)) {
    SPLICING_ERROR("Invalid matrix-vector sizes", SPLICING_EINVAL);
  }
  if (splicing_vector_size(y) != (transpose ? n : m)) {
    SPLICING_ERROR("Invalid matrix-vector sizes", SPLICING_EINVAL);
  }    

  splicingdgemv_(&trans, &m, &n, &alpha, VECTOR(a->data), &m,
		 VECTOR(*x), &inc, &beta, VECTOR(*y), &inc);
  
  return 0;
}
