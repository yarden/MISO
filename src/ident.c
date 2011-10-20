
#include "splicing_error.h"
#include "splicing.h"

/* Might mess up 'exstart', 'exend' and 'exidx'!!! */

int splicing_numeric_cigar(splicing_vector_int_t *exstart, 
			   splicing_vector_int_t *exend,
			   splicing_vector_int_t *exidx,
			   int noiso, size_t genestart, 
			   splicing_vector_int_t *result, 
			   int skip) {
			   
  int i;

  genestart += skip;

  splicing_vector_int_clear(result);
  SPLICING_CHECK(splicing_vector_int_push_back(result, 0));
  for (i=0; i<noiso; i++) {
    int pos=VECTOR(*exidx)[i];
    int pos2=VECTOR(*exidx)[i+1];
    
    /* We don't start from the beginning of the gene */
    if (VECTOR(*exstart)[pos] < genestart) { 
      while (VECTOR(*exstart)[pos] < genestart) {
	if (VECTOR(*exend)[pos] < genestart) { 
	  /* skip whole exon */
	  pos++;
	} else {
	  /* skip first part of exon */
	  VECTOR(*exstart)[pos] = genestart;
	}
      }
    }
    
    /* Does the isoform start where we start? */
    if (pos < pos2 && VECTOR(*exstart)[pos] > genestart) { 
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

