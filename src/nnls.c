
#include "splicing.h"
#include "splicing_error.h"

int splicing_nnls(splicing_matrix_t *A, splicing_vector_t *B, 
		  splicing_vector_t *X, double *rnorm, 
		  splicing_vector_long_t *index, long int *nsetp) {
  
  long int m=splicing_matrix_nrow(A);
  long int n=splicing_matrix_ncol(A);
  splicing_vector_t w, zz;
  long int mode;

  if (splicing_vector_size(B) != m) { 
    SPLICING_ERROR("Invalid matrix dimensions", SPLICING_EINVAL);
  }

  SPLICING_CHECK(splicing_vector_init(&w, n));
  SPLICING_FINALLY(splicing_vector_destroy, &w);
  SPLICING_CHECK(splicing_vector_init(&zz, m));
  SPLICING_FINALLY(splicing_vector_destroy, &zz);
  SPLICING_CHECK(splicing_vector_resize(X, n));
  SPLICING_CHECK(splicing_vector_long_resize(index, n));

  nnls_(&MATRIX(*A, 0, 0), &m, &m, &n, VECTOR(*B), VECTOR(*X), rnorm,
	VECTOR(w), VECTOR(zz), VECTOR(*index), &mode, nsetp);
  
  splicing_vector_destroy(&zz);
  splicing_vector_destroy(&w);
  SPLICING_FINALLY_CLEAN(2);
  
  if (mode==1) {
    return 0;
  } else {
    SPLICING_ERROR("NNLS did not converge", SPLICING_DIVERGED);
  }
}

