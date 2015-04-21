
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "splicing_matrix.h"
#include "splicing_error.h"
#include "splicing.h"
#include "splicing_random.h"

#define SIGMA (0.2/noiso/noiso)

int splicing_drift_proposal_init_paired(int noiso, int noChains, 
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
				 double normalVar, double numDevs) {

  SPLICING_CHECK(splicing_matrix_resize(respsi, noiso, noChains));
  SPLICING_CHECK(splicing_matrix_resize(resalpha, noiso-1, noChains));

  switch (start) {
  case SPLICING_MISO_START_AUTO:
    if (noiso != 2) {
      int i, j;
      for (i=0; i<noiso; i++) {	
	for (j=0; j<noChains; j++) {
	  MATRIX(*respsi, i, j) = 1.0/noiso; 
	}
      }
      for (i=0; i<noiso-1; i++) { 
	for (j=0; j<noChains; j++) {
	  MATRIX(*resalpha, i, j) = 1.0/(noiso-1);
	}
      }
      *ressigma = SIGMA;
    } else {
      int j;
      for (j=0; j<noChains; j++) {
	MATRIX(*respsi, 0, j) = RNG_UNIF01();
	MATRIX(*respsi, 1, j) = 1 - MATRIX(*respsi, j, 0);
	MATRIX(*resalpha, 0, j) = 0.0;
      }
      *ressigma = SIGMA;
    }
    break;
  case SPLICING_MISO_START_UNIFORM:
    {
      int i, j;
      for (i=0; i<noiso; i++) {	
	for (j=0; j<noChains; j++) {
	  MATRIX(*respsi, i, j) = 1.0/noiso; 
	}
      }
      for (i=0; i<noiso-1; i++) { 
	for (j=0; j<noChains; j++) {
	  MATRIX(*resalpha, i, j) = 0;
	}
      }
      *ressigma = SIGMA;
    }
    break;
  case SPLICING_MISO_START_RANDOM:
    {
      splicing_vector_t alpha, tmp;
      int i, j;
      SPLICING_CHECK(splicing_vector_init(&alpha, noiso));
      SPLICING_FINALLY(splicing_vector_destroy, &alpha);
      for (i=0; i<noiso; i++) { VECTOR(alpha)[i] = 1.0; }
      for (j=0; j<noChains; j++) {
	splicing_vector_view(&tmp, &MATRIX(*respsi, 0, j), noiso);
	SPLICING_CHECK(splicing_rng_get_dirichlet(&splicing_rng_default,
						  &alpha, &tmp));
      }
      SPLICING_CHECK(splicing_logit(respsi, resalpha, noiso-1, noChains));
      splicing_vector_destroy(&alpha);
      SPLICING_FINALLY_CLEAN(1);
    }
    break;
  case SPLICING_MISO_START_GIVEN:
    splicing_matrix_update(respsi, start_psi);
    SPLICING_CHECK(splicing_logit(respsi, resalpha, noiso-1, noChains));
    *ressigma = SIGMA;
    break;
  case SPLICING_MISO_START_LINEAR:
    {
      splicing_vector_t tmp, tmp2;
      int j;
      SPLICING_CHECK(splicing_vector_init(&tmp, noiso));
      SPLICING_FINALLY(splicing_vector_destroy, &tmp);
      if (!paired) { 
	SPLICING_CHECK(splicing_solve_gene(gff, gene, readLength, 
				       overHang, position, cigarstr,
				       /*match_matrix=*/ 0, /*nomatch=*/ 0,
				       /*assignment_matrix=*/ 0,
				       &tmp, /*residuals=*/ 0,
				       /*scale=*/ 1));
      } else {
	SPLICING_CHECK(splicing_solve_gene_paired(gff, gene, readLength,
				       overHang, position,
				       cigarstr, fragmentProb, fragmentStart,
				       normalMean, normalVar, numDevs,
				       /*match_matrix=*/ 0, /*nomatch=*/ 0,
				       /*assignment_matrix=*/ 0, &tmp,
				       /*residuals=*/ 0, /*scale*/ 1));
      }
      
      for (j=0; j<noChains; j++) {
	splicing_vector_view(&tmp2, &MATRIX(*respsi, 0, j), noiso);
	splicing_vector_update(&tmp2, &tmp);
      }
      SPLICING_CHECK(splicing_logit(respsi, resalpha, noiso-1,
				    noChains));
      splicing_vector_destroy(&tmp);
      SPLICING_FINALLY_CLEAN(1);
      *ressigma = SIGMA;
    }
    break;
  }
  
  return 0;
}

int splicing_drift_proposal_propose_paired(int noiso, int noChains, 
				    const splicing_matrix_t *alpha,
				    double sigma, 
				    splicing_matrix_t *respsi,
				    splicing_matrix_t *resalpha) {

  int len=noiso-1;
  double sumpsi=0.0;
  int i;
  
  SPLICING_CHECK(splicing_matrix_resize(respsi, len+1, noChains));
  SPLICING_CHECK(splicing_mvrnorm(alpha, sigma, resalpha, len));
  SPLICING_CHECK(splicing_logit_inv(resalpha, respsi, len, noChains));
  for (i=0; i<noChains; i++) {
    splicing_vector_t col;
    splicing_vector_view(&col, &MATRIX(*respsi, 0, i), len);
    sumpsi = splicing_vector_sum(&col);
    /* SPLICING_CHECK(splicing_vector_resize(respsi, len+1)); */
    MATRIX(*respsi, len, i) = 1-sumpsi;
  }
  
  return 0;
 }
