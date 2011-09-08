
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
int R_splicing_SEXP_to_matrix_int(SEXP pm, splicing_matrix_int_t *m);
int R_splicing_SEXP_to_gff(SEXP pgff, splicing_gff_t *gff);
int R_splicing_SEXP_to_strvector(SEXP pv, splicing_strvector_t *v);
int R_splicing_SEXP_to_exons(SEXP pexons, splicing_vector_int_t *exons);
int R_splicing_SEXP_to_isoforms(SEXP piso, splicing_vector_int_t *iso);

/* These are not needed in Python, SAM/BAM is handled differently there */

typedef struct splicing_reads_t {
  int noPairs, noSingles, paired;
  splicing_strvector_t chrname;
  splicing_vector_int_t chrlen;
  splicing_vector_int_t chr;
  splicing_strvector_t qname;
  splicing_strvector_t cigar;
  splicing_vector_int_t position;
  splicing_vector_int_t flags;
  splicing_vector_int_t pairpos;
  splicing_vector_int_t mapq;
  splicing_vector_int_t rnext;
  splicing_vector_int_t tlen;
  splicing_strvector_t seq;
  splicing_strvector_t qual;
  splicing_vector_int_t mypair;
} splicing_reads_t;

typedef enum { SPLICING_SAMBAM_AUTO, 
	       SPLICING_SAMBAM_SAM,
	       SPLICING_SAMBAM_BAM } splicing_sambam_type_t;

int splicing_reads_init(splicing_reads_t *reads);
void splicing_reads_destroy(splicing_reads_t *reads);
int splicing_read_sambam(const char *filename,
			 splicing_sambam_type_t filetype,
			 splicing_reads_t *reads);
int splicing_read_sambam_region(const char *filename,
				const char *indexfile,
				splicing_sambam_type_t filetype,
				const char *region,
				splicing_reads_t *reads);

SEXP R_splicing_reads_to_SEXP(const splicing_reads_t *reads);

#endif
