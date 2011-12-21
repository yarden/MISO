
#ifndef SPLICING_H
#define SPLICING_H

#include <stdio.h>
#include <stdlib.h>

#include "splicing_vector.h"
#include "splicing_matrix.h"
#include "splicing_lapack.h"
#include "splicing_random.h"

double round(double);
double fmin(double, double);

void splicing_qsort(void *base, size_t nel, size_t width,
		    int (*compar)(const void *, const void *));
void splicing_qsort_r(void *base, size_t nel, size_t width, void *thunk,
		      int (*compar)(void *, const void *, const void *));

#define SPLICING_NA_INTEGER (-1)
#define SPLICING_NA_REAL    (-1.0)

typedef struct {
  size_t size;
  size_t asize;
  char **table;
  int free;
} splicing_strvector_t;

extern const char *splicing_strvector_zero;
#define SPLICING_STRVECTOR_ZERO ((char*) splicing_strvector_zero)

int splicing_strvector_init(splicing_strvector_t *v, size_t size);
void splicing_strvector_destroy(splicing_strvector_t *v);
size_t splicing_strvector_size(const splicing_strvector_t *v);
int splicing_strvector_append(splicing_strvector_t *v, const char *str);
int splicing_strvector_append2(splicing_strvector_t *v, const char *str, 
			       size_t len);
int splicing_strvector_reserve(splicing_strvector_t *v, size_t size);
const char *splicing_strvector_get(const splicing_strvector_t *v, 
				   size_t idx);

int splicing_strvector_search(const splicing_strvector_t *v, 
			      const char *key, size_t *idx);
int splicing_strvector_fprint(const splicing_strvector_t *v, FILE *file);
int splicing_strvector_print(const splicing_strvector_t *v);
int splicing_strvector_clear(splicing_strvector_t *v);
int splicing_strvector_permute(splicing_strvector_t *v, 
			       const splicing_vector_t *idx);
int splicing_strvector_ipermute(splicing_strvector_t *v, 
				const splicing_vector_int_t *idx);

double splicing_dnorm(double x, double mu, double sigma);
double splicing_logdnorm(double x, double mu, double sigma);

extern const char *splicing_types[];

typedef enum { SPLICING_TYPE_GENE, SPLICING_TYPE_MRNA, 
	       SPLICING_TYPE_EXON, SPLICING_TYPE_CDS, 
	       SPLICING_TYPE_START_CODON, SPLICING_TYPE_STOP_CODON }
  splicing_type_t;

typedef enum { SPLICING_STRAND_PLUS, SPLICING_STRAND_MINUS, 
	       SPLICING_STRAND_UNKNOWN } splicing_strand_t;

/* TODO: arbitrary attributes */

typedef struct {
  size_t n;
  splicing_strvector_t seqids;
  splicing_strvector_t sources;

  splicing_vector_int_t genes;
  splicing_vector_int_t transcripts;
  splicing_vector_int_t seqid;
  splicing_vector_int_t source;

  splicing_vector_int_t type;
  splicing_vector_int_t start;
  splicing_vector_int_t end;
  splicing_vector_t score;
  splicing_vector_int_t strand;
  splicing_vector_int_t phase;
  splicing_strvector_t ID;
  splicing_vector_int_t parent;

  int nogenes, notranscripts;

  const char *last_gene_id, *last_mrna_id, *last_seqid, *last_source;
  int last_gene_no, last_mrna_no;  
} splicing_gff_t;

int splicing_gff_init(splicing_gff_t *gff, size_t size);
void splicing_gff_destroy(splicing_gff_t *gff);
void splicing_gff_destroy2(void *gff);
int splicing_gff_reserve(splicing_gff_t *gff, size_t size);
size_t splicing_gff_size(const splicing_gff_t *gff);
int splicing_gff_append(splicing_gff_t *gff, const char *seqid, 
			const char *source, splicing_type_t type, int start,
			int end, double score, splicing_strand_t strand, 
			int phase, const char *ID, const char *parent);
int splicing_gff_read(FILE *input, splicing_gff_t *gff);
int splicing_gff_write(FILE *output, const splicing_gff_t *gff);
int splicing_gff_noiso_one(const splicing_gff_t *gff, size_t gene, 
			   size_t *noiso);
int splicing_gff_isolength_one(const splicing_gff_t *gff, size_t gene, 
			       splicing_vector_int_t *isolength);
int splicing_gff_noiso(const splicing_gff_t *gff, 
		       splicing_vector_int_t *noiso);
int splicing_gff_isolength(const splicing_gff_t *gff,
			   splicing_vector_int_t *isolength,
			   splicing_vector_int_t *isolength_idx);
int splicing_gff_nogenes(const splicing_gff_t *gff, size_t *nogenes);

int splicing_gff_noexons_one(const splicing_gff_t *gff, size_t gene,
			     splicing_vector_int_t *noexons);

int splicing_gff_exon_start_end(const splicing_gff_t *gff, 
				splicing_vector_int_t *start,
				splicing_vector_int_t *end,
				splicing_vector_int_t *idx, int gene);

int splicing_gff_gene_start_end_one(const splicing_gff_t *gff, size_t gene,
				    size_t *start, size_t *end);
int splicing_gff_gene_start_end(const splicing_gff_t *gff, 
				splicing_vector_int_t *start,
				splicing_vector_int_t *end);

int splicing_gff_fprint_gene(const splicing_gff_t *gff, 
			     FILE *outfile, int gene);
int splicing_gff_print_gene(const splicing_gff_t *gff, 
			    int gene);
int splicing_gff_fprint(const splicing_gff_t *gff, 
			FILE *outfile);
int splicing_gff_print(const splicing_gff_t *gff);

typedef struct splicing_miso_rundata_t {
  int noIso, noIters, maxIters, noBurnIn, noLag, noAccepted, noRejected,
    noChains, noSamples;
} splicing_miso_rundata_t;

typedef enum splicing_miso_start_t {
  SPLICING_MISO_START_AUTO,
  SPLICING_MISO_START_UNIFORM,
  SPLICING_MISO_START_RANDOM,
  SPLICING_MISO_START_GIVEN,
  SPLICING_MISO_START_LINEAR } splicing_miso_start_t;

typedef enum splicing_miso_stop_t {
  SPLICING_MISO_STOP_FIXEDNO,
  SPLICING_MISO_STOP_CONVERGENT_MEAN,
} splicing_miso_stop_t;

int splicing_matchIso(const splicing_gff_t *gff, int gene, 
		      const splicing_vector_int_t *position,
		      const char **cigarstr, int overHang, int readLength,
		      splicing_matrix_t *result);

int splicing_matchIso_paired(const splicing_gff_t *gff, int gene,
			     const splicing_vector_int_t *position,
			     const char **cigarstr, int readLength,
			     int overHang, 
			     const splicing_vector_t *fragmentProb,
			     int fragmentStart, double normalMean,
			     double normalVar, double numDevs,
			     splicing_matrix_t *result,
			     splicing_matrix_int_t *fragmentLengths);

int splicing_parse_cigar(const char **cigar, size_t noreads,
			 splicing_vector_int_t *numcigar,
			 splicing_vector_int_t *cigaridx, 
			 splicing_vector_int_t *cigarlength,
			 int maxReadLength);

int splicing_order_matches(const splicing_matrix_t *matches,
			   splicing_vector_int_t *order);

int splicing_i_miso_classes(const splicing_matrix_t *match_matrix,
			    const splicing_vector_int_t *match_order,
			    splicing_matrix_t *class_templates,
			    splicing_vector_t *class_counts, 
			    splicing_matrix_t *bin_class_templates,
			    splicing_vector_t *bin_class_counts);

int splicing_i_check_convergent_mean(splicing_matrix_t *chainMeans, 
				     splicing_matrix_t *chainVars, 
				     const splicing_matrix_t *samples,
				     int *shouldstop);

int splicing_miso(const splicing_gff_t *gff, size_t gene,
		  const splicing_vector_int_t *position,
		  const char **cigarstr, int readLength, int overHang,
		  int noChains, int noIterations, int maxIterations, 
		  int noBurnIn, int noLag, const splicing_vector_t *hyperp, 
		  splicing_miso_start_t start, splicing_miso_stop_t stop,
		  const splicing_matrix_t *start_psi,
		  splicing_matrix_t *samples, splicing_vector_t *logLik, 
		  splicing_matrix_t *match_matrix,
		  splicing_matrix_t *class_templates,
		  splicing_vector_t *class_counts,
		  splicing_vector_int_t *assignment,
		  splicing_miso_rundata_t *rundata);

int splicing_miso_paired(const splicing_gff_t *gff, size_t gene,
			 const splicing_vector_int_t *position,
			 const char **cigarstr, int readLength, int overHang,
			 int noChains, int noIterations, int maxIterations,
			 int noBurnIn, int noLag, 
			 const splicing_vector_t *hyperp, 
			 splicing_miso_start_t start, 
			 splicing_miso_stop_t stop, 
			 const splicing_matrix_t *start_psi,
			 const splicing_vector_t *fragmentProb, 
			 int fragmentStart,
			 double normalMean, double normalVar, double numDevs,
			 splicing_matrix_t *samples, 
			 splicing_vector_t *logLik, 
			 splicing_matrix_t *match_matrix,
			 splicing_matrix_t *class_templates,
			 splicing_vector_t *class_counts,
			 splicing_matrix_t *bin_class_templates,
			 splicing_vector_t *bin_class_count,
			 splicing_vector_int_t *assignment,
			 splicing_miso_rundata_t *rundata);

int splicing_reassign_samples(const splicing_matrix_t *matches, 
			      const splicing_vector_int_t *match_order,
			      const splicing_matrix_t *psi, 
			      int noiso, int noChains, 
			      splicing_matrix_int_t *result);

int splicing_mvplogisnorm(const splicing_vector_t *theta, 
			  const splicing_vector_t *mu, 
			  double sigma, int len, double *score);

int splicing_score_iso(const splicing_vector_t *psi, int noiso, 
		       const splicing_vector_int_t *assignment, int noreads,
		       const splicing_vector_int_t *peffisolen, double *res);

int splicing_ldirichlet(const splicing_vector_t *x, 
			const splicing_vector_t *alpha, int len, 
			double *res);

int splicing_mvrnorm(const splicing_matrix_t *mu, double sigma, 
		     splicing_matrix_t *resalpha, int len);

int splicing_logit_inv(const splicing_matrix_t *x, 
		       splicing_matrix_t *res, int len, int noChains);

int splicing_score_joint(const splicing_matrix_int_t *assignment,
			 int no_reads, int noChains, 
			 const splicing_matrix_t *psi, 
			 const splicing_vector_t *hyper, 
			 const splicing_vector_int_t *effisolen,
			 const splicing_vector_t *isoscores, 
			 splicing_vector_t *score);

int splicing_drift_proposal_init(int noiso, int noChains, 
				 splicing_matrix_t *respsi, 
				 splicing_matrix_t *resalpha,
				 double *ressigma,
				 splicing_miso_start_t start,
				 const splicing_matrix_t *start_psi,
				 const splicing_gff_t *gff, int gene, 
				 int readLength, int overHang, 
				 const splicing_vector_int_t *position,
				 const char **cigarstr, int paired, 
				 const splicing_vector_t *fragmentProb,
				 int fragmentStart, double normalMean,
				 double normalVar, double numDevs);

int splicing_drift_proposal_propose(int noiso, int noChains, 
				    const splicing_matrix_t *alpha,
				    double sigma, 
				    splicing_matrix_t *respsi,
				    splicing_matrix_t *resalpha);

int splicing_drift_proposal_score(int noiso, int noChains, 
				  const splicing_matrix_t *psi,
				  const splicing_matrix_t *otheralpha, 
				  double sigma,
				  splicing_vector_t *resscore);

int splicing_metropolis_hastings_ratio(const splicing_matrix_int_t *ass,
				       int no_reads, int noChains,
				       const splicing_matrix_t *psiNew,
				       const splicing_matrix_t *alphaNew,
				       const splicing_matrix_t *psi, 
				       const splicing_matrix_t *alpha,
				       double sigma,
				       int noiso, 
				       const splicing_vector_int_t *effisolen,
				       const splicing_vector_t *hyperp, 
				       const splicing_vector_t *isoscores,
				       int full, splicing_vector_t *acceptP, 
				       splicing_vector_t *pcJS, 
				       splicing_vector_t *ppJS);

int splicing_assignment_matrix(const splicing_gff_t *gff, size_t gene,
			       int readLength, int overHang, 
			       splicing_matrix_t *matrix);

int splicing_numeric_cigar(const splicing_vector_int_t *exstart, 
			   const splicing_vector_int_t *exend,
			   const splicing_vector_int_t *exidx,
			   int noiso, size_t genestart, 
			   splicing_vector_int_t *result);

int splicing_create_gene(const splicing_vector_int_t *exons,
			 const splicing_vector_int_t *isoforms,
			 const char *id, const char *seqid, 
			 const char *source, splicing_strand_t strand, 
			 splicing_gff_t *extend);

int splicing_simulate_reads(const splicing_gff_t *gff, int gene,
			    const splicing_vector_t *expression,
			    int noreads, int readLength,
			    splicing_vector_int_t *isoform, 
			    splicing_vector_int_t *position, 
			    splicing_strvector_t *cigar, 
			    splicing_vector_t *sample_prob);

int splicing_normal_fragment(double normalMean, double normalVar, 
			     double numDevs, int minLength,
			     splicing_vector_t *fragmentProb,
			     int *fragmentStart);

int splicing_simulate_paired_reads(const splicing_gff_t *gff, int gene,
				   const splicing_vector_t *expression,
				   int noreads, int readLength,
				   const splicing_vector_t *fragmentProb,
				   int fragmentStart, double normalMean,
				   double normalVar, double numDevs,
				   splicing_vector_int_t *isoform,
				   splicing_vector_int_t *position,
				   splicing_strvector_t *cigar, 
				   splicing_vector_t *sampleprob);

int nnls_(double *a, long int *mda, long int *m, long int *n, double *b,
	  double *x, double *rnorm, double *w, double *zz, long int *index, 
	  long int *mode, long int *nsetp);

int splicing_nnls(splicing_matrix_t *A, splicing_vector_t *B, 
		  splicing_vector_t *X, double *rnorm, 
		  splicing_vector_long_t *index, long int *nsetp);

int splicing_solve_gene(const splicing_gff_t *gff, size_t gene, 
			int readLength, int overHang,
			const splicing_vector_int_t *position, 
			const char **cigarstr,
			splicing_matrix_t *match_matrix,
			splicing_matrix_t *assignment_matrix, 
			splicing_vector_t *expression);

int splicing_solve_gene_paired(const splicing_gff_t *gff, size_t gene,
			       int readLength, int overHang,
			       const splicing_vector_int_t *position,
			       const char **cigarstr,
			       const splicing_vector_t *fragmentProb,
			       int fragmentStart, double normalMean,
			       double normalVar, double numDevs,
			       splicing_matrix_t *match_matrix,
			       splicing_matrix_t *assignment_matrix,
			       splicing_vector_t *expression);

typedef struct splicing_gff_converter_t {
  size_t noiso;
  splicing_vector_int_t exstart, exend, exidx, exlim, shift;
} splicing_gff_converter_t;

int splicing_gff_converter_init(const splicing_gff_t *gff, size_t gene,
				splicing_gff_converter_t *converter);

void splicing_gff_converter_destroy(splicing_gff_converter_t *converter);

int splicing_iso_to_genomic(const splicing_gff_t *gff, size_t gene, 
			    const splicing_vector_int_t *isoform,
			    const splicing_gff_converter_t *converter,
			    splicing_vector_int_t *position);

int splicing_iso_to_genomic_1(const splicing_gff_t *gff, size_t gene,
			      int isoform, int position, 
			      const splicing_gff_converter_t *converter,
			      int *result);

int splicing_iso_to_genomic_all(const splicing_gff_t *gff, size_t gene,
				int position, 
				const splicing_gff_converter_t *converter,
				splicing_vector_int_t *result);

int splicing_genomic_to_iso(const splicing_gff_t *gff, size_t gene,
			    const splicing_vector_int_t *position, 
			    const splicing_gff_converter_t *converter,
			    splicing_matrix_int_t *isopos);

int splicing_genomic_to_iso_1(const splicing_gff_t *gff, size_t gene,
			      int isoform, int position, 
			      const splicing_gff_converter_t *converter,
			      int *result);

int splicing_genomic_to_iso_all(const splicing_gff_t *gff, size_t gene,
				int position, 
				const splicing_gff_converter_t *converter,
				splicing_vector_int_t *result);

int splicing_paired_assignment_matrix(const splicing_gff_t *gff, size_t gene,
				      int readLength, int overHang,
				      const splicing_vector_t *fragmentProb,
				      int fragmentStart, double normalMean,
				      double normalVar, double numDevs,
				      splicing_matrix_t *matrix);

int splicing_paired_assignment_matrix_old(const splicing_gff_t *gff, 
				  size_t gene, int readLength, int overHang,
				  const splicing_vector_t *fragmentProb,
				  int fragmentStart, double normalMean,
				  double normalVar, double numDevs,
				  splicing_matrix_t *matrix);

int splicing_reassign_samples_paired(
			     const splicing_matrix_t *matches, 
			     const splicing_vector_int_t *match_order,
			     const splicing_matrix_t *psi, 
			     int noiso, int noChains, int fragmentStart, 
			     splicing_matrix_int_t *result);

int splicing_score_iso_paired(const splicing_vector_t *psi, int noiso, 
			      const splicing_vector_int_t *assignment, 
			      const splicing_vector_int_t *pisolen, 
			      const splicing_vector_t *assscores,
			      double *res);

int splicing_score_joint_paired(const splicing_matrix_int_t *assignment,
				int no_reads, int noChains,
				const splicing_matrix_t *psi, 
				const splicing_vector_t *hyper, 
				const splicing_vector_int_t *isolen,
				const splicing_matrix_t *isoscores, 
				const splicing_vector_t *assscores,
				const splicing_matrix_int_t *fragmentLength,
				int fragmentStart, 
				splicing_vector_t *score);

int splicing_metropolis_hastings_ratio_paired(
			      const splicing_matrix_int_t *ass,
			      int no_reads, int noChains,
			      const splicing_matrix_t *psiNew,
			      const splicing_matrix_t *alphaNew,
			      const splicing_matrix_t *psi, 
			      const splicing_matrix_t *alpha,
			      double sigma, int noiso, 
			      const splicing_vector_int_t *isolen,
			      const splicing_vector_t *hyperp, 
			      const splicing_matrix_t *isoscores,
			      const splicing_vector_t *assscores,
			      const splicing_matrix_int_t *fragmentLength,
			      int fragmentStart, int full, 
			      splicing_vector_t *acceptP, 
			      splicing_vector_t *pcJS, 
			      splicing_vector_t *ppJS);

int splicing_dgesdd(const splicing_matrix_t *matrix, 
		    splicing_vector_t *values);

typedef enum { SPLICING_COMPLEXITY_RELATIVE,
	       SPLICING_COMPLEXITY_ABSOLUTE } splicing_complexity_t;

typedef enum { SPLICING_NORM_2, 
	       SPLICING_NORM_1, 
	       SPLICING_NORM_INFINITY } splicing_norm_t;

int splicing_gene_complexity(const splicing_gff_t *gff, size_t gene,
			     int readLength, int overHang, 
			     splicing_complexity_t type,
			     splicing_norm_t norm, int paired, 
			     const splicing_vector_t *fragmentProb,
			     int fragmentStart, double normalMean, 
			     double normalVar, double numDevs,
			     double *complexity);

int splicing_rng_get_dirichlet(splicing_rng_t *rng, 
			       const splicing_vector_t *alpha, 
			       splicing_vector_t *result);

#endif
