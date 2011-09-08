
#include "splicing.h"
#include "splicing_error.h"

int splicingdgesdd_(char *jobz, int *m, int *n, double *a,
		    int *lda, double *s, double *u, int *ldu, 
		    double *vt, int *ldvt, double *work, int *lwork, 
		    int *iwork, int *info);

int splicing_dgesdd(const splicing_matrix_t *matrix, 
		    splicing_vector_t *values) {

  splicing_matrix_t tmp;
  int m=splicing_matrix_nrow(matrix);
  int n=splicing_matrix_ncol(matrix);
  int lda=m, minmn= m < n ? m : n, maxmn = m < n ? n : m;
  int lwork=-1;
  int info=0;
  splicing_vector_t work;
  splicing_vector_int_t iwork;
  char jobz='N';
  int dummy=1;
  double dummy2;
  
  SPLICING_CHECK(splicing_matrix_copy(&tmp, matrix));
  SPLICING_FINALLY(splicing_matrix_destroy, &tmp);
  SPLICING_CHECK(splicing_vector_init(&work, 1));
  SPLICING_FINALLY(splicing_vector_destroy, &work);
  SPLICING_CHECK(splicing_vector_int_init(&iwork, 8*minmn));
  SPLICING_FINALLY(splicing_vector_int_destroy, &iwork);

  SPLICING_CHECK(splicing_vector_resize(values, minmn));

  /* Get the optiomal lwork first*/
  splicingdgesdd_(&jobz, &m, &n, &MATRIX(tmp,0,0), &lda, VECTOR(*values),
		  /*U=*/ &dummy2, /*LDU=*/ &dummy, 
		  /*VT=*/ &dummy2, /*LDVT=*/ &dummy, 
		  VECTOR(work), &lwork, VECTOR(iwork), &info);

  lwork = VECTOR(work)[0];
  SPLICING_CHECK(splicing_vector_resize(&work, lwork));

  /* Now do the SVD */
  splicingdgesdd_(&jobz, &m, &n, &MATRIX(tmp,0,0), &lda, VECTOR(*values),
		  /*U=*/ &dummy2, /*LDU=*/ &dummy, 
		  /*VT=*/ &dummy2, /*LDVT=*/ &dummy, 
		  VECTOR(work), &lwork, VECTOR(iwork), &info);

  if (info != 0) { 
    SPLICING_ERROR("Cannot calculate SVD", SPLICING_ELAPACK);
  }

  splicing_vector_destroy(&work);
  splicing_vector_int_destroy(&iwork);
  splicing_matrix_destroy(&tmp);
  SPLICING_FINALLY_CLEAN(3);
  
  return 0;
}
