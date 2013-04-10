
#include "splicing_error.h"
#include "splicing.h"

/* Might mess up 'exstart', 'exend' and 'exidx'!!! */

int splicing_numeric_cigar(const splicing_vector_int_t *exstart, 
			   const splicing_vector_int_t *exend,
			   const splicing_vector_int_t *exidx,
			   int noiso, size_t genestart, 
			   splicing_vector_int_t *result) {
			   
  int i;

  splicing_vector_int_clear(result);
  SPLICING_CHECK(splicing_vector_int_push_back(result, 0));
  for (i=0; i<noiso; i++) {
    int pos=VECTOR(*exidx)[i];
    int pos2=VECTOR(*exidx)[i+1];
    
    /* Does the isoform start where we start? */
    if (VECTOR(*exstart)[pos] != genestart) { 
      SPLICING_CHECK(splicing_vector_int_push_back(result, 
					genestart-VECTOR(*exstart)[pos]));
    }
    while (pos<pos2) {
      int l=VECTOR(*exend)[pos] - VECTOR(*exstart)[pos] + 1;
      SPLICING_CHECK(splicing_vector_int_push_back(result, l));
      if (pos+1 < pos2) {
	int l=VECTOR(*exend)[pos]-VECTOR(*exstart)[pos+1] + 1;
	SPLICING_CHECK(splicing_vector_int_push_back(result, l));
      }
      pos++;
    }
    SPLICING_CHECK(splicing_vector_int_push_back(result, 0));
  }
  
  return 0;
}

