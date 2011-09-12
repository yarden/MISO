
#include "splicing.h"
#include "splicing_error.h"

int splicing_gene_complexity(const splicing_gff_t *gff, size_t gene,
			     int readLength, int overHang,
			     splicing_complexity_t type,
			     splicing_norm_t norm, int paired,
			     const splicing_vector_t *fragmentProb,
			     int fragmentStart, double normalMean, 
			     double normalVar, double numDevs,
			     double *complexity) {
  
  splicing_matrix_t assignment_matrix;

  SPLICING_CHECK(splicing_matrix_init(&assignment_matrix, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &assignment_matrix);

  if (!paired) {
    SPLICING_CHECK(splicing_assignment_matrix(gff, gene, readLength, 
					      overHang, &assignment_matrix));
  } else {
    SPLICING_CHECK(splicing_paired_assignment_matrix(gff, gene, readLength, 
						     overHang, fragmentProb, 
						     fragmentStart,
						     normalMean, normalVar,
						     numDevs, 
						     &assignment_matrix));
  }

  switch (type) {
  case SPLICING_COMPLEXITY_RELATIVE:
    switch (norm) {
      splicing_vector_t values;
      int i, n;

    case SPLICING_NORM_2:

      SPLICING_CHECK(splicing_vector_init(&values, 0));
      SPLICING_FINALLY(splicing_vector_destroy, &values);
      SPLICING_CHECK(splicing_dgesdd(&assignment_matrix, &values));
      n=splicing_vector_size(&values);
      for (i=n-1; i>=0 && VECTOR(values)[i] < 1e-14; i--) ;
      *complexity = VECTOR(values)[0] / VECTOR(values)[i];
      splicing_vector_destroy(&values);
      SPLICING_FINALLY_CLEAN(1);
      break;

    case SPLICING_NORM_1:

      SPLICING_ERROR("One norm not implemented", SPLICING_UNIMPLEMENTED);
      break;

    case SPLICING_NORM_INFINITY:

      SPLICING_ERROR("Infinity norm not implemented", SPLICING_UNIMPLEMENTED);
      break;

    }
    break;
  case SPLICING_COMPLEXITY_ABSOLUTE:
    SPLICING_ERROR("Absolute complexity not implemented", 
		   SPLICING_UNIMPLEMENTED);
    break;
  }

  splicing_matrix_destroy(&assignment_matrix);
  SPLICING_FINALLY_CLEAN(1);

  return 0;
}
