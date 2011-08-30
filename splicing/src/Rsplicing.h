
#ifndef R_SPLICING_H
#define R_SPLICING_H

#include <R.h>
#include <Rdefines.h>

#include "splicing.h"

SEXP R_splicing_getListElement(SEXP list, const char *str);

SEXP R_splicing_vector_int_to_SEXP(const splicing_vector_int_t *v);
SEXP R_splicing_vector_to_SEXP(const splicing_vector_t *v);
SEXP R_splicing_matrix_to_SEXP(const splicing_matrix_t *m);
SEXP R_splicing_matrix_int_to_SEXP(const splicing_matrix_int_t *m);
SEXP R_splicing_gff_to_SEXP(splicing_gff_t *gff);
SEXP R_splicing_strvector_to_SEXP(const splicing_strvector_t *strv);
SEXP R_splicing_vector_int_index_to_SEXP(const splicing_vector_int_t *v,
					 const splicing_vector_int_t *idx);
SEXP R_splicing_miso_rundata_to_SEXP(const splicing_miso_rundata_t *data);

int R_splicing_SEXP_to_vector(SEXP pv, splicing_vector_t *v);
int R_splicing_SEXP_to_vector_int(SEXP pv, splicing_vector_int_t *v);
int R_splicing_SEXP_to_matrix(SEXP pm, splicing_matrix_t *m);
int R_splicing_SEXP_to_gff(SEXP pgff, splicing_gff_t *gff);
int R_splicing_SEXP_to_strvector(SEXP pv, splicing_strvector_t *v);
int R_splicing_SEXP_to_exons(SEXP pexons, splicing_vector_int_t *exons);
int R_splicing_SEXP_to_isoforms(SEXP piso, splicing_vector_int_t *iso);

int splicing_read_sambam(const char *filename, 
			 splicing_strvector_t *chrname, 
			 splicing_vector_int_t *chrlen,
			 splicing_vector_int_t *chr,
			 splicing_strvector_t *qname,
			 splicing_strvector_t *cigar,
			 splicing_vector_int_t *position, 
			 splicing_vector_int_t *flag,
			 splicing_vector_int_t *pairpos,
			 int *noPairs, int *noSingles, int *paired,
			 splicing_vector_int_t *mapq,
			 splicing_vector_int_t *rnext,
			 splicing_vector_int_t *tlen,
			 splicing_strvector_t *seq,
			 splicing_strvector_t *qual);


#endif
