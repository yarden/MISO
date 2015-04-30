
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "splicing_matrix.h"
#include "splicing_error.h"
#include "splicing.h"
#include "splicing_random.h"
#include "splicing_memory.h"

#define CUMSUM() do {							\
    int j;								\
    double *cptr = curr;						\
    for (j=0, noValid=0, sumpsi=0.0; j<noiso; j++, cptr++) {		\
      if (*cptr != 0) {							\
	sumpsi += MATRIX(*psi, j, k);					\
	VECTOR(validIso)[noValid] = j;					\
	VECTOR(cumsum)[noValid] = sumpsi;				\
	noValid++;							\
      }									\
    }									\
  } while (0)

int splicing_reassign_samples1(const splicing_matrix_t *matches,
			       const splicing_vector_int_t *match_order,
			       const splicing_matrix_t *psi, 
			       int noiso, int noChains, 
			       splicing_matrix_int_t *result);

int splicing_reassign_samples(const splicing_vector_ptr_t *all_matches, /* matrix_t */
			      const splicing_vector_ptr_t *all_match_order, /* vector_int_t */
			      const splicing_vector_ptr_t *all_psi, /* matrix_t */
			      int noiso, int noChains, 
			      splicing_vector_ptr_t *all_result) { /* matrix_int_t */

  int rep, norep = (int) splicing_vector_ptr_size(all_matches);
  for (rep = 0; rep < norep; rep++) {
    const splicing_matrix_t *matches = VECTOR(*all_matches)[rep];
    const splicing_vector_int_t *match_order = VECTOR(*all_match_order)[rep];
    const splicing_matrix_t *psi = VECTOR(*all_psi)[rep];
    splicing_matrix_int_t *result = VECTOR(*all_result)[rep];
    splicing_reassign_samples1(matches, match_order, psi, noiso, noChains,
			       result);
  }

  return 0;
}  

  

/* TODO: we could actually speed this up, by storing the 
   indices where the cumulative sum vector must be updated. 
   These indices do not change at all. With this we can spare
   the memcmp operations. */

int splicing_reassign_samples1(const splicing_matrix_t *matches,
			       const splicing_vector_int_t *match_order,
			       const splicing_matrix_t *psi, 
			       int noiso, int noChains, 
			       splicing_matrix_int_t *result) {
    
  int noreads = (int) splicing_matrix_ncol(matches);
  int i, k, w;
  double *prev, *curr;
  double rand, sumpsi;
  int noValid;
  int *order=VECTOR(*match_order);
  splicing_vector_t cumsum;
  splicing_vector_int_t validIso;  

  SPLICING_CHECK(splicing_vector_init(&cumsum, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &cumsum);
  SPLICING_CHECK(splicing_vector_int_init(&validIso, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &validIso);

  SPLICING_CHECK(splicing_matrix_int_resize(result, noreads, noChains));

  if (noreads == 0) { return 0; }  

  for (k=0; k<noChains; k++) {

    prev = curr = &MATRIX(*matches, 0, order[0]);
    CUMSUM();

    for (i=0; i<noreads; i++) {
      curr = &MATRIX(*matches, 0, order[i]);

      /* Maybe we need to update the cumulative sum */
      if (memcmp(prev, curr, sizeof(double)*noiso) != 0) { CUMSUM(); }
      
      if (noValid == 0) {
	MATRIX(*result, order[i], k) = -1;
      } else if (noValid == 1) {
	MATRIX(*result, order[i], k) = VECTOR(validIso)[0];
      } else if (noValid == 2) { 
	rand = RNG_UNIF01() * sumpsi;
	w = (rand < VECTOR(cumsum)[0]) ? VECTOR(validIso)[0] : 
	  VECTOR(validIso)[1];
	MATRIX(*result, order[i], k) = w;
      } else {
	/* Draw */
	rand = RNG_UNIF01() * sumpsi;
	/* TODO: Binary search for interval, if many classes */
	for (w=0; rand > VECTOR(cumsum)[w]; w++) ;
	MATRIX(*result, order[i], k) = VECTOR(validIso)[w];
      }

      prev=curr;
    }
  }

  splicing_vector_int_destroy(&validIso);
  splicing_vector_destroy(&cumsum);
  SPLICING_FINALLY_CLEAN(2);

  return 0;
}

/* We only handle a special case here, where sigma is a diagonal
   matrix with identical elements. In this case it is easy to 
   invert it, or calculate its determinant. */

int splicing_mvplogisnorm(const splicing_vector_t *theta, 
			  const splicing_vector_t *mu, 
			  double sigma, int len, double *score) {
  
  double covarConst = pow(2 * M_PI * sigma, -0.5 * len);
  double ltheta=1.0, prodTheta=1.0, expPart=0.0, pdfVal;
  int i;
  
  for (i=0; i<len; i++) {
    double at=VECTOR(*theta)[i];
    ltheta -= at;
    prodTheta *= at;
  }
  prodTheta = 1.0 / prodTheta / ltheta;
  
  for (i=0; i<len; i++) {
    double tmp=log(VECTOR(*theta)[i] / ltheta) - VECTOR(*mu)[i];
    expPart += (-0.5) * tmp * tmp / sigma;
  }
  
  pdfVal = covarConst * prodTheta * exp(expPart);

  *score = log(pdfVal);
  
  return 0;
}

int splicing_score_iso(const splicing_vector_t *psi, int noiso, 
		       const splicing_vector_int_t *assignment, int noreads,
		       const splicing_vector_int_t *peffisolen, double *res) {
  int *effisolen = VECTOR(*peffisolen);
  double sum, maxpsieff, score;
  splicing_vector_t logpsi;
  int i;

  SPLICING_CHECK(splicing_vector_init(&logpsi, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &logpsi);

  /* Calculate the normalization factor */
  VECTOR(logpsi)[0] = log(VECTOR(*psi)[0]) + log(effisolen[0]);
  for (maxpsieff=VECTOR(logpsi)[0], i=1; i<noiso; i++) {
    VECTOR(logpsi)[i] = log(VECTOR(*psi)[i]) + log(effisolen[i]);
    if (VECTOR(logpsi)[i] > maxpsieff) { maxpsieff = VECTOR(logpsi)[i]; }
  }
  for (sum=0.0, i=0; i<noiso; i++) {
    sum += exp(VECTOR(logpsi)[i]-maxpsieff);
  }
  sum = log(sum) + maxpsieff;
  
  /* Normalize */
  for (i=0; i<noiso; i++) {
    VECTOR(logpsi)[i] -= sum;
  }
  
  /* Calculate score, based on assignments */
  for (score=0.0, i=0; i<noreads; i++) {
    if (VECTOR(*assignment)[i] != -1) {
      score += VECTOR(logpsi)[ VECTOR(*assignment)[i] ];
    }
  }

  splicing_vector_destroy(&logpsi);
  SPLICING_FINALLY_CLEAN(1);

  *res = score;
  return 0;
}

int splicing_ldirichlet(const splicing_vector_t *x, 
			const splicing_vector_t *alpha, int len, 
			double *res) {
  double score=0.0, alphasum=0.0, lgalphasum=0.0;
  int i;

  for (i=0; i<len; i++) {
    alphasum += VECTOR(*alpha)[i];
    lgalphasum += lgamma(VECTOR(*alpha)[i]);
    score += (VECTOR(*alpha)[i] - 1.0) * log(VECTOR(*x)[i]);
  }
  
  score += lgamma(alphasum);
  score -= lgalphasum;

  *res = score;
  return 0;
}

double splicing_logit1(double v) {
  return - log( (1 - v) / v );
}

int splicing_llogistic(const splicing_vector_t *x,
		       double logistic_mean, double logistic_var,
		       double *res) {

  /* We know that we have two isoforms only, and only deal with the
     first isoform, the second is just 1 - x[0] */
  double logit_psi = splicing_logit1(VECTOR(*x)[0]);

  *res = splicing_logdnorm(logit_psi, logistic_mean, sqrt(logistic_var));

  return 0;
}

int splicing_mvrnorm(const splicing_matrix_t *mu, double sigma, 
		     splicing_matrix_t *resalpha, int len) {
  int i, j;
  int noChains = (int) splicing_matrix_ncol(mu);
  double sqrtsigma = len == 1 ? sigma : sqrt(sigma);

  SPLICING_CHECK(splicing_matrix_resize(resalpha, len, noChains));

  for (j=0; j<noChains; j++) {
    for (i=0; i<len; i++) {
      MATRIX(*resalpha, i, j) = MATRIX(*mu, i, j) + 
	sqrtsigma * RNG_NORMAL(0,1);
    }
  }

  return 0;
}

int splicing_logit(const splicing_matrix_t *x,
		   splicing_matrix_t *res, int len, int noChains) {

  int i, j;
    
  SPLICING_CHECK(splicing_matrix_resize(res, len, noChains));
  
  for (j=0; j<noChains; j++) {
    double logyn=log(MATRIX(*x, len, j));
    for (i=0; i<len; i++) {
      MATRIX(*res, i, j) = log(MATRIX(*x, i, j)) - logyn;
    }
  }
  
  return 0;
}

int splicing_logit_inv(const splicing_matrix_t *x, 
		       splicing_matrix_t *res, int len, int noChains) {
  int i, j;

  if (splicing_matrix_nrow(res) <  len || 
      splicing_matrix_ncol(res) != noChains) { 
    SPLICING_ERROR("`res' has an illegal size", SPLICING_EINVAL);
  }

  for (j=0; j<noChains; j++) {
    double sumexp=0.0;
    for (i=0; i<len; i++) {
      sumexp += exp(MATRIX(*x, i, j));
    }
    sumexp += 1.0;
  
    for (i=0; i<len; i++) {
      MATRIX(*res, i, j) = exp(MATRIX(*x, i, j)) / sumexp;
    }
  }
  
  return 0;
}

int splicing_score_joint(splicing_algorithm_t algorithm,
			 const splicing_matrix_int_t *assignment,
			 int no_reads, int noChains, 
			 const splicing_matrix_t *psi, 
			 splicing_miso_hyperprior_t *hyperprior,
			 const splicing_vector_int_t *effisolen,
			 const splicing_vector_t *isoscores,
			 const splicing_matrix_t *match,
			 const splicing_matrix_t *assignmentMatrix,
			 const splicing_vector_t *matches,
			 splicing_vector_t *score) {

  int i, j, noiso = (int) splicing_vector_int_size(effisolen);

  SPLICING_CHECK(splicing_vector_resize(score, noChains));

  for (j=0; j<noChains; j++) {
    double readProb = 0.0, assProb = 0.0, psiProb;
    splicing_vector_t tmp;
    splicing_vector_int_t tmp2;
    splicing_vector_view(&tmp, &MATRIX(*psi, 0, j), noiso);

    if (algorithm == SPLICING_ALGO_REASSIGN) {
      /* Scores the reads */
      for (i=0; i<no_reads; i++) {
	if (MATRIX(*assignment, i, j) != -1) {
	  readProb += VECTOR(*isoscores)[ MATRIX(*assignment, i, j) ];
	}
      }
    } else if (algorithm == SPLICING_ALGO_MARGINAL) {
      /* TODO: we could replace this with a matrix-vector multiplication */
      for (i=0; i<no_reads; i++) {
	double isoscore=0.0;
	int k;
	for (k=0; k<noiso; k++) {
	  isoscore += MATRIX(*match, k, i) * MATRIX(*psi, k, j);
	}
	if (isoscore != 0) { readProb += log(isoscore); }
      }
    } else if (algorithm == SPLICING_ALGO_CLASSES) {
      int no_classes = (int) splicing_matrix_ncol(assignmentMatrix);
      /* TODO: we could replace this with a matrix-vector multiplication */
      for (i=0; i<no_classes; i++) {
	double score=0.0;
	int k;
	for (k=0; k<noiso; k++) {
	  score += MATRIX(*assignmentMatrix, k, i) * MATRIX(*psi, k, j);
	}
	if (score != 0) { readProb += log(score) * VECTOR(*matches)[i]; }
      }
    }
    
    /* Score isoforms */
    if (algorithm == SPLICING_ALGO_REASSIGN) {
      splicing_vector_int_view(&tmp2, &MATRIX(*assignment, 0, j), no_reads);
      SPLICING_CHECK(splicing_score_iso(&tmp, noiso, &tmp2, no_reads,
					effisolen, &assProb));
    }
    if (hyperprior->prior == SPLICING_MISO_PRIOR_DIRICHLET) {
      SPLICING_CHECK(splicing_ldirichlet(
         &tmp, &(hyperprior->dirichlet_hyperp), noiso, &psiProb));
    } else if (hyperprior->prior == SPLICING_MISO_PRIOR_LOGISTIC) {
      SPLICING_CHECK(splicing_llogistic(
	 &tmp, hyperprior->logistic_mean, hyperprior->logistic_var,
	 &psiProb));
    } else {
      SPLICING_ERROR("Unknown prior", SPLICING_EINVAL);
    }

    VECTOR(*score)[j] = readProb + assProb + psiProb;
  }

  return 0;
}

int splicing_rng_get_dirichlet(splicing_rng_t *rng, 
			       const splicing_vector_t *alpha, 
			       splicing_vector_t *result) {

  int i, l = (int) splicing_vector_size(alpha);
  double sum=0.0;
  
  SPLICING_CHECK(splicing_vector_resize(result, l));
  for (i=0; i<l; i++) { 
    VECTOR(*result)[i] = RNG_GAMMA(VECTOR(*alpha)[i], 1.0);
    sum += VECTOR(*result)[i];
  }
  for (i=0; i<l; i++) {
    VECTOR(*result)[i] /= sum;
  }

  return 0;
}

int splicing_drift_proposal_init1(int noiso, int noChains, 
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

int splicing_drift_proposal_init(int noiso, int noChains, 
				 splicing_vector_ptr_t *all_psi,
				 splicing_vector_ptr_t *all_alpha,
				 double *ressigma,
				 splicing_miso_start_t start,
				 const splicing_matrix_t *start_psi,
				 const splicing_gff_t *gff, int gene, 
				 int readLength, int overHang,
				 const splicing_replicate_reads_t *reads,
				 int paired, 
				 const splicing_vector_t *fragmentProb,
				 int fragmentStart, double normalMean,
				 double normalVar, double numDevs) {

  int i, noReplicates = (int) splicing_vector_ptr_size(all_psi);

  for (i = 0; i < noReplicates; i++) {
    splicing_matrix_t *psi = VECTOR(*all_psi)[i];
    splicing_matrix_t *alpha = VECTOR(*all_alpha)[i];
    const splicing_vector_int_t *position =
      splicing_replicate_reads_pos(reads, i);
    const char **cigarstr = splicing_replicate_reads_cigar(reads, i);
    splicing_drift_proposal_init1(noiso, noChains, psi, alpha, ressigma,
				  start, start_psi, gff, gene, readLength,
				  overHang, position, cigarstr, paired,
				  fragmentProb, fragmentStart, normalMean,
				  normalVar, numDevs);
  }
       
  return 0;
}


#define SIGMA (0.2/noiso/noiso)

int splicing_drift_proposal_init1(int noiso, int noChains, 
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

int splicing_drift_proposal_propose1(int noiso, int noChains, 
				     const splicing_matrix_t *alpha,
				     double sigma, 
				     splicing_matrix_t *respsi,
				     splicing_matrix_t *resalpha);

int splicing_drift_proposal_propose(int noiso, int noChains, 
				    const splicing_vector_ptr_t *all_alpha,
				    double sigma, 
				    splicing_vector_ptr_t *all_respsi, /* matrix_t */
				    splicing_vector_ptr_t *all_resalpha) { /* matrix_t */

  int i, noReplicates = (int) splicing_vector_ptr_size(all_respsi);
  for (i = 0; i < noReplicates; i++) {
    const splicing_matrix_t *alpha = VECTOR(*all_alpha)[i];
    splicing_matrix_t *respsi = VECTOR(*all_respsi)[i];
    splicing_matrix_t *resalpha = VECTOR(*all_resalpha)[i];
    splicing_drift_proposal_propose1(noiso, noChains, alpha, sigma,
				     respsi, resalpha);
  }

  return 0;
}

int splicing_drift_proposal_propose1(int noiso, int noChains, 
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

int splicing_drift_proposal_score(int noiso, int noChains, 
				  const splicing_matrix_t *psi,
				  const splicing_matrix_t *otheralpha, 
				  double sigma,
				  splicing_vector_t *resscore) {

  double score;
  int i;
  SPLICING_CHECK(splicing_vector_resize(resscore, noChains));
  for (i=0; i<noChains; i++) {
    splicing_vector_t col1, col2;
    splicing_vector_view(&col1, &MATRIX(*psi, 0, i), noiso);
    splicing_vector_view(&col2, &MATRIX(*otheralpha, 0, i), noiso-1);
    SPLICING_CHECK(splicing_mvplogisnorm(&col1, &col2, sigma, noiso-1, 
					 &score));
    VECTOR(*resscore)[i]=score;
  }
  return 0;
}

int splicing_metropolis_hastings_ratio(splicing_algorithm_t algorithm,
				       const splicing_matrix_int_t *ass,
				       int no_reads, int noChains,
				       const splicing_matrix_t *psiNew,
				       const splicing_matrix_t *alphaNew,
				       const splicing_matrix_t *psi, 
				       const splicing_matrix_t *alpha,
				       double sigma,
				       int noiso, 
				       const splicing_matrix_t *match,
				       const splicing_matrix_t *assignment,
				       const splicing_vector_t *matches,
				       const splicing_vector_int_t *effisolen,
				       splicing_miso_hyperprior_t *hyperprior,
				       const splicing_vector_t *isoscores,
				       int full, splicing_vector_t *acceptP, 
				       splicing_vector_t *pcJS, 
				       splicing_vector_t *ppJS) {

  int i;
  splicing_vector_t ptoCS, ctoPS;

  SPLICING_CHECK(splicing_vector_resize(acceptP, noChains));
  SPLICING_CHECK(splicing_vector_resize(pcJS, noChains));
  SPLICING_CHECK(splicing_vector_resize(ppJS, noChains));

  SPLICING_CHECK(splicing_vector_init(&ptoCS, noChains));
  SPLICING_FINALLY(splicing_vector_destroy, &ptoCS);
  SPLICING_CHECK(splicing_vector_init(&ctoPS, noChains));
  SPLICING_FINALLY(splicing_vector_destroy, &ctoPS);

  SPLICING_CHECK(splicing_score_joint(algorithm, ass, no_reads, noChains,
				      psiNew, hyperprior,
				      effisolen, isoscores, match, assignment,
				      matches, ppJS));
  SPLICING_CHECK(splicing_score_joint(algorithm, ass, no_reads, noChains, psi,
				      hyperprior, effisolen, isoscores, match,
				      assignment, matches, pcJS));
  
  SPLICING_CHECK(splicing_drift_proposal_score(noiso, noChains, psi, alphaNew,
					       sigma, &ptoCS));
  SPLICING_CHECK(splicing_drift_proposal_score(noiso, noChains, psiNew, alpha,
					       sigma, &ctoPS));
  
  if (full) {
    for (i=0; i<noChains; i++) {
      VECTOR(*acceptP)[i] = exp(VECTOR(*ppJS)[i] + VECTOR(ptoCS)[i] - 
				(VECTOR(*pcJS)[i] + VECTOR(ctoPS)[i]));
    }
  } else {
    for (i=0; i<noChains; i++) {
      VECTOR(*acceptP)[i] = exp(VECTOR(*ppJS)[i] - VECTOR(*pcJS)[i]);
    }
  }

  splicing_vector_destroy(&ctoPS);
  splicing_vector_destroy(&ptoCS);
  SPLICING_FINALLY_CLEAN(2);

  return 0;
}

/* chainMeans and chainVars must have the correct size! */

int splicing_i_check_convergent_mean(splicing_matrix_t *chainMeans, 
				     splicing_matrix_t *chainVars, 
				     const splicing_matrix_t *samples,
				     int *shouldstop) {

  /* Check convergence. See Gelman, Carlin, Stern and Rubin:
     Bayesian Data Analysis, pp. 296, 2nd edition, for
     details. */
  int i /* sample */, j /* chain */, k /* isoform */,
    l /* sample per chain */;
  splicing_vector_t B, W, mean, rhat;
  int noiso = (int) splicing_matrix_nrow(chainMeans);
  int noChains = (int) splicing_matrix_ncol(chainMeans);
  int noSamples = (int) splicing_matrix_ncol(samples);

  splicing_matrix_null(chainVars);
  memcpy(&MATRIX(*chainMeans, 0, 0), &MATRIX(*samples, 0, 0), 
	 noChains * noiso * sizeof(double));
  
  for (i=noChains, j=0, l=1; i<noSamples; i++, j = (j+1) % noChains) {
    for (k=0; k<noiso; k++) {
      double mk=MATRIX(*chainMeans, k, j) + 
	(MATRIX(*samples, k, i) - MATRIX(*chainMeans, k, j)) / l;
      double sk=MATRIX(*chainVars, k, j) + 
	(MATRIX(*samples, k, i) - MATRIX(*chainMeans, k, j)) *
	(MATRIX(*samples, k, i) - mk);
      MATRIX(*chainMeans, k, j) = mk;
      MATRIX(*chainVars, k, j) = sk;
    }
    if (j==noChains-1) { l++; }
  }
  
  SPLICING_CHECK(splicing_vector_init(&B, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &B);
  SPLICING_CHECK(splicing_vector_init(&W, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &W);
  SPLICING_CHECK(splicing_vector_init(&mean, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &mean);
  SPLICING_CHECK(splicing_vector_init(&rhat, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &rhat);
  
  for (k=0; k<noiso; k++) { 
    for (j=0; j<noChains; j++) { 
      VECTOR(mean)[k] += MATRIX(*chainMeans, k, j);
    }
    VECTOR(mean)[k] /= noChains;
  }
  
  for (k=0; k<noiso; k++) {
    for (j=0; j<noChains; j++) { 
      double t=MATRIX(*chainMeans, k, j) - VECTOR(mean)[k];
      VECTOR(B)[k] += t*t;
    }
    VECTOR(B)[k] *= noSamples / (noChains-1.0);
  }
  
  for (k=0; k<noiso; k++) {
    for (j=0; j<noChains; j++) {
      double t=MATRIX(*chainVars, k, j);
      VECTOR(W)[k] += t * t;
    }
    VECTOR(W)[k] /= noChains;
  }
  
  for (k=0; k<noiso; k++) {
    VECTOR(rhat)[k] = 
      sqrt(((noSamples - 1.0)/noSamples * VECTOR(W)[k] + VECTOR(B)[k] /
	    noSamples) / VECTOR(W)[k]);
  }
  
  for (k=0, *shouldstop=1; k<noiso; k++) {
    *shouldstop = *shouldstop && VECTOR(rhat)[k] <= 1.1;
  }
  splicing_vector_destroy(&rhat);
  splicing_vector_destroy(&mean);
  splicing_vector_destroy(&W);
  splicing_vector_destroy(&B);
  SPLICING_FINALLY_CLEAN(4);
  
  return 0;
}

/* Replicates: we need to go over all replicate individually,
   for the most part. Replicates are only considered together
   when drawing population mean and variance. */

int splicing_miso(
		  /* INPUT */

		  const splicing_gff_t *gff, size_t gene,
		  const splicing_replicate_reads_t *reads,
		  int readLength, int overHang,
		  int noChains, int noIterations, 
		  int maxIterations, int noBurnIn, int noLag,
		  splicing_miso_hyperprior_t *hyperprior,
		  splicing_algorithm_t algorithm,
		  splicing_miso_start_t start,
		  splicing_miso_stop_t stop,
		  const splicing_matrix_t *start_psi,

		  /* OUTPUT */

		  /* The samples for the psi vectors, for each replicate */
		  splicing_vector_ptr_t *samples,	  /* matrix_t */

		  /* Log-likelihood for the samples */
		  splicing_vector_t *logLik,

		  /* Which isoforms the reads match, for each replicate */
		  splicing_vector_ptr_t *match_matrix,    /* matrix_t */

		  /* Which isoforms the various classes are matching */
		  splicing_matrix_t *class_templates,

		  /* Number of reads in various classes, per replicate */
		  splicing_vector_ptr_t *class_counts,

		  /* An example assignment of reads to isoforms */
		  splicing_vector_ptr_t *assignment,	  /* vector_int_t */

		  /* Information about the MCMC */
		  splicing_miso_rundata_t *rundata) {

  /* These are gene dependent, and the same for all replicates */

  size_t noiso;
  int noSamples = noChains * (noIterations - noBurnIn) / noLag;  
  int i, j, m=0, lagCounter=0, noS=0;
  splicing_vector_int_t effisolen;
  splicing_vector_t isoscores;
  splicing_vector_int_t noexons;
  int shouldstop=0;
  splicing_matrix_t assignmentMatrix;
  int noClasses;
  splicing_vector_int_t noReads;
  int noReplicates = splicing_replicate_reads_noreps(reads);
  
  /* These are reused for each replicate */

  splicing_vector_t acceptP, cJS, pJS;
  double sigma;
  
  /* These are different for each replicate */

  splicing_vector_ptr_t *mymatch_matrix = match_matrix,
    vmatch_matrix;			       /* matrix_t */
  splicing_vector_ptr_t match_order;	       /* vector_int_t */
  splicing_vector_ptr_t chainMeans, chainVars; /* matrix_t */
  splicing_vector_ptr_t matches;	       /* vector_t */
  splicing_vector_ptr_t vass;		       /* matrix_int_t, 
						  assignments to isoforms */
  int have_vass = 0;

  /* All these contain matrix_t */
  splicing_vector_ptr_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;

  if (algorithm != SPLICING_ALGO_REASSIGN &&
      algorithm != SPLICING_ALGO_MARGINAL &&
      algorithm != SPLICING_ALGO_CLASSES) {
    SPLICING_ERROR("`algorithm` is invalid", SPLICING_EINVAL);
  }

  if (start == SPLICING_MISO_START_GIVEN && !start_psi) {
    SPLICING_ERROR("`start_psi' must be given when "
		   "starting from a given PSI", SPLICING_EINVAL);
  }

  if ( (class_templates ? 1 : 0) + (class_counts ? 1 : 0) == 1) {
    SPLICING_ERROR("Only one of `class_templates' and `class_counts' is "
		   "given", SPLICING_EINVAL);
  }

  if (overHang==0) { overHang=1; }
  if (overHang < 1 || overHang >= readLength / 2) {
    SPLICING_ERROR("Overhang length invalid. Must be between 0 and "
		   "readLength/2", SPLICING_EINVAL);
  }

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));

  if (hyperprior->prior != SPLICING_MISO_PRIOR_AUTO &&
      hyperprior->prior != SPLICING_MISO_PRIOR_DIRICHLET &&
      hyperprior->prior != SPLICING_MISO_PRIOR_LOGISTIC) {
    SPLICING_ERROR("'prior' is invalid", SPLICING_EINVAL);
  }

  if (hyperprior->prior == SPLICING_MISO_PRIOR_AUTO) {
    if (noiso == 2) {
      hyperprior->prior = SPLICING_MISO_PRIOR_LOGISTIC;
    } else {
      hyperprior->prior = SPLICING_MISO_PRIOR_DIRICHLET;
    }
  }

  if (hyperprior->prior == SPLICING_MISO_PRIOR_DIRICHLET &&
      splicing_vector_size(&(hyperprior->dirichlet_hyperp)) != noiso) {
    SPLICING_ERROR("Invalid hyperparameter vector length", 
		   SPLICING_EINVAL);
  }

  if (hyperprior->prior == SPLICING_MISO_PRIOR_LOGISTIC && noiso != 2) {
    SPLICING_ERROR("Logistic prior currently only works for two isoforms",
		   SPLICING_UNIMPLEMENTED);
  }

  if (hyperprior->logistic_var < 0) {
    SPLICING_ERROR("Variance of the logistic prior must be non-negative",
		   SPLICING_EINVAL);
  }

  if (noChains < 1) { 
    SPLICING_ERROR("Number of chains must be at least one.", 
		   SPLICING_EINVAL);
  }

  if (stop==SPLICING_MISO_STOP_CONVERGENT_MEAN && noChains == 1) {
    SPLICING_ERROR("Cannot access convergence with one chain only", 
		   SPLICING_EINVAL);
  }

  if (start_psi && 
      (splicing_matrix_nrow(start_psi) != noiso ||
       splicing_matrix_ncol(start_psi) != noChains)) {
    SPLICING_ERROR("Given PSI has wrong size", SPLICING_EINVAL);
  }

  rundata->noIso=(int) noiso;
  rundata->noIters=noIterations;
  rundata->noBurnIn=noBurnIn;
  rundata->noLag=noLag;
  rundata->noAccepted = rundata->noRejected = 0;
  rundata->noChains = noChains;
  rundata->noSamples = noSamples;

  SPLICING_CHECK(splicing_vector_init(&acceptP, noChains));
  SPLICING_FINALLY(splicing_vector_destroy, &acceptP);
  SPLICING_CHECK(splicing_vector_init(&cJS, noChains));
  SPLICING_FINALLY(splicing_vector_destroy, &cJS);
  SPLICING_CHECK(splicing_vector_init(&pJS, noChains));
  SPLICING_FINALLY(splicing_vector_destroy, &pJS);

  SPLICING_CHECK(splicing_vector_int_init(&noReads, noReplicates));
  for (i = 0; i < noReplicates; i++) {
    VECTOR(noReads)[i] = splicing_replicate_reads_noreads(reads, i);
  }

  /* We need the assignment for the classic algorithm, or if they
     were requested. */
  if (algorithm == SPLICING_ALGO_REASSIGN || assignment) {
    have_vass = 1;
    SPLICING_CHECK(splicing_vector_ptr_init(&vass, noReplicates));
    splicing_vector_ptr_set_item_destructor(&vass,
	(splicing_finally_func_t *) splicing_matrix_int_destroy_free);
    SPLICING_FINALLY(splicing_vector_ptr_destroy, &vass);
    for (i = 0; i < noReplicates; i++) {
      splicing_matrix_int_t *vass1 =
	splicing_Calloc(1, splicing_matrix_int_t);
      if (!vass1) {
	SPLICING_ERROR("No memory for assignment", SPLICING_ENOMEM);
      }
      VECTOR(vass)[i] = vass1;
      SPLICING_CHECK(splicing_matrix_int_init(vass1, VECTOR(noReads)[i],
					      noChains));
    }
  }

  SPLICING_CHECK(splicing_vector_ptr_init(&vpsi, noReplicates));
  splicing_vector_ptr_set_item_destructor(&vpsi,
    (splicing_finally_func_t *) splicing_matrix_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &vpsi);

  SPLICING_CHECK(splicing_vector_ptr_init(&vpsiNew, noReplicates));
  splicing_vector_ptr_set_item_destructor(&vpsiNew,
    (splicing_finally_func_t *) splicing_matrix_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &vpsiNew);

  SPLICING_CHECK(splicing_vector_ptr_init(&valpha, noReplicates));
  splicing_vector_ptr_set_item_destructor(&valpha,
    (splicing_finally_func_t *) splicing_matrix_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &valpha);

  SPLICING_CHECK(splicing_vector_ptr_init(&valphaNew, noReplicates));
  splicing_vector_ptr_set_item_destructor(&valphaNew,
    (splicing_finally_func_t *) splicing_matrix_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &valphaNew);

  for (i = 0; i < noReplicates; i++) {
    splicing_matrix_t *psi = splicing_Calloc(1, splicing_matrix_t);
    splicing_matrix_t *psinew = splicing_Calloc(1, splicing_matrix_t);
    splicing_matrix_t *alpha = splicing_Calloc(1, splicing_matrix_t);
    splicing_matrix_t *alphanew = splicing_Calloc(1, splicing_matrix_t);
    if (!psi || !psinew || !alpha || !alphanew) {
      SPLICING_ERROR("No memory to run MISO", SPLICING_ENOMEM);
    }
    VECTOR(vpsi)[i] = psi;
    VECTOR(vpsiNew)[i] = psinew;
    VECTOR(valpha)[i] = alpha;
    VECTOR(valphaNew)[i] = alphanew;
    SPLICING_CHECK(splicing_matrix_init(psi, noiso, noChains));
    SPLICING_CHECK(splicing_matrix_init(psinew, noiso, noChains));
    SPLICING_CHECK(splicing_matrix_init(alpha, noiso-1, noChains));
    SPLICING_CHECK(splicing_matrix_init(alphanew, noiso-1, noChains));
  }
  
  if (match_matrix) {
    for (i = 0; i < noReplicates; i++) {
      splicing_matrix_t *mat = VECTOR(*match_matrix)[i];
      SPLICING_CHECK(splicing_matrix_resize(mat, noiso, VECTOR(noReads)[i]));
    }
  } else {
    mymatch_matrix = &vmatch_matrix;
    SPLICING_CHECK(splicing_vector_ptr_init(mymatch_matrix, noReplicates));
    splicing_vector_ptr_set_item_destructor(mymatch_matrix,
      (splicing_finally_func_t *) splicing_matrix_destroy_free);
    SPLICING_FINALLY(splicing_vector_ptr_destroy, mymatch_matrix);
    for (i = 0; i < noReplicates; i++) {
      splicing_matrix_t *mat = splicing_Calloc(1, splicing_matrix_t);
      if (!mat) SPLICING_ERROR("No memory for match matrix", SPLICING_ENOMEM);
      VECTOR(*mymatch_matrix)[i] = mat;
      SPLICING_CHECK(splicing_matrix_init(mat, noiso, VECTOR(noReads)[i]));
    }
  }
  SPLICING_CHECK(splicing_vector_ptr_init(&match_order, noReplicates));
  splicing_vector_ptr_set_item_destructor(&match_order,
    (splicing_finally_func_t *) splicing_vector_int_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &match_order);
  for (i = 0; i < noReplicates; i++) {
    splicing_vector_int_t *vec = splicing_Calloc(1, splicing_vector_int_t);
    if (!vec) { SPLICING_ERROR("No memory for match order", SPLICING_ENOMEM); }
    VECTOR(match_order)[i] = vec;
    SPLICING_CHECK(splicing_vector_int_init(vec, VECTOR(noReads)[i]));
  }

  for (i = 0; i < noReplicates; i++) {
    const splicing_vector_int_t *position =
      splicing_replicate_reads_pos(reads, i);
    const char **cigarstr = splicing_replicate_reads_cigar(reads, i);
    splicing_matrix_t *mm = VECTOR(*mymatch_matrix)[i];
    splicing_vector_int_t *mo = VECTOR(match_order)[i];
    SPLICING_CHECK(splicing_matchIso(gff, (int) gene, position, cigarstr,
				     overHang, readLength, mm));
    SPLICING_CHECK(splicing_order_matches(mm, mo));
  }
  
  if (class_templates && class_counts) {
    SPLICING_CHECK(splicing_vector_ptr_resize(class_counts, noReplicates));
    splicing_vector_ptr_set_item_destructor(class_counts,
      (splicing_finally_func_t *) splicing_vector_destroy_free);
    SPLICING_FINALLY(splicing_vector_ptr_destroy, class_counts);
    for (i = 0; i < noReplicates; i++) {
      splicing_matrix_t *mm = VECTOR(*mymatch_matrix)[i];
      splicing_vector_int_t *mo = VECTOR(match_order)[i];
      splicing_vector_t *cc = splicing_Calloc(1, splicing_vector_t);
      if (!cc) { SPLICING_ERROR("No memory for class counts", SPLICING_ENOMEM); }
      VECTOR(*class_counts)[i] = cc;
      SPLICING_CHECK(splicing_i_miso_classes(mm, mo, class_templates, cc,
					     /*bin_class_templates=*/ 0,
					     /*bin_class_counts=*/ 0));
    }
  }

  SPLICING_CHECK(splicing_vector_int_init(&effisolen, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &effisolen);
  SPLICING_CHECK(splicing_vector_init(&isoscores, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &isoscores);
  SPLICING_CHECK(splicing_gff_isolength_one(gff, gene, &effisolen));
  SPLICING_CHECK(splicing_vector_int_init(&noexons, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &noexons);
  SPLICING_CHECK(splicing_gff_noexons_one(gff, gene, &noexons));
  for (i=0; i<noiso; i++) { 
    int nox=VECTOR(noexons)[i];
    /* The following is only approximate if there are some short exons
       and overHang is not one */
    int l=VECTOR(effisolen)[i] - readLength+1 - 2*(nox-1)*(overHang-1);
    VECTOR(effisolen)[i] = l > 0 ? l : 0;
    VECTOR(isoscores)[i] = -log((double) l);
  }
  splicing_vector_int_destroy(&noexons);
  SPLICING_FINALLY_CLEAN(1);

  /* Calculate assignment matrix and the number of reads in each class,
     if we are running that algorithm */
  if (algorithm == SPLICING_ALGO_CLASSES) {
    SPLICING_CHECK(splicing_matrix_init(&assignmentMatrix, 0, 0));
    SPLICING_FINALLY(splicing_matrix_destroy, &assignmentMatrix);
    SPLICING_CHECK(splicing_assignment_matrix(gff, gene, readLength, overHang,
					      &assignmentMatrix));
    splicing_matrix_norm_row(&assignmentMatrix);
    noClasses = (int) splicing_matrix_ncol(&assignmentMatrix);
    SPLICING_CHECK(splicing_vector_ptr_init(&matches, noReplicates));
    splicing_vector_ptr_set_item_destructor(&matches,
      (splicing_finally_func_t*) splicing_vector_destroy_free);
    SPLICING_FINALLY(splicing_vector_ptr_destroy, &matches);
    for (i = 0; i < noReplicates; i++) {
      const splicing_vector_int_t *position =
	splicing_replicate_reads_pos(reads, i);
      const char **cigarstr = splicing_replicate_reads_cigar(reads, i);
      splicing_matrix_t *mm = VECTOR(*mymatch_matrix)[i];
      splicing_vector_t *matches1 = splicing_Calloc(1, splicing_vector_t);
      if (!matches1) { SPLICING_ERROR("No memory for matches", SPLICING_ENOMEM); }
      VECTOR(matches)[i] = matches1;
      SPLICING_CHECK(splicing_vector_init(matches1, noClasses));
      SPLICING_CHECK(splicing_getMatchVector(gff, (int) gene,
					     VECTOR(noReads)[i], position,
					     cigarstr, overHang, readLength,
					     mm, &assignmentMatrix,
					     matches1));
    }
  }

  /* For the marginal algorithm, the match_matrix will contain
     probabilities divided by effective isoform length */
  if (algorithm == SPLICING_ALGO_MARGINAL) {
    int r;
    for (r=0; i<noReplicates; r++) {
      splicing_matrix_t *mm = VECTOR(*mymatch_matrix)[r];
      int nor = VECTOR(noReads)[r];
      for (i=0; i<noiso; i++) {
	for (j=0; j<nor; j++) {
	  if (VECTOR(effisolen)[i] != 0) {
	    MATRIX(*mm, i, j) /= VECTOR(effisolen)[i];
	  }
	}
      }
    }
  }

  SPLICING_CHECK(splicing_vector_ptr_init(&chainMeans, noReplicates));
  splicing_vector_ptr_set_item_destructor(&chainMeans,
    (splicing_finally_func_t *) splicing_matrix_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &chainMeans);
  for (i = 0; i < noReplicates; i++) {
    splicing_matrix_t *cm = splicing_Calloc(1, splicing_matrix_t);
    if (!cm) { SPLICING_ERROR("No memory for chain means", SPLICING_ENOMEM); }
    VECTOR(chainMeans)[i] = cm;
    SPLICING_CHECK(splicing_matrix_init(cm, noiso, noChains));
  }

  SPLICING_CHECK(splicing_vector_ptr_init(&chainVars, noReplicates));
  splicing_vector_ptr_set_item_destructor(&chainVars,
    (splicing_finally_func_t *) splicing_matrix_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &chainVars);
  for (i = 0; i < noReplicates; i++) {
    splicing_matrix_t *cv = splicing_Calloc(1, splicing_matrix_t);
    if (!cv) { SPLICING_ERROR("No memory for chain vars", SPLICING_ENOMEM); }
    VECTOR(chainVars)[i] = cv;
    SPLICING_CHECK(splicing_matrix_init(cv, noiso, noChains));
  }

  SPLICING_CHECK(splicing_vector_ptr_resize(samples, noReplicates));
  splicing_vector_ptr_set_item_destructor(samples,
    (splicing_finally_func_t *) splicing_matrix_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, samples);
  for (i = 0; i < noReplicates; i++) {
    splicing_matrix_t *ss = splicing_Calloc(1, splicing_matrix_t);
    if (!ss) { SPLICING_ERROR("No memory for samples", SPLICING_ENOMEM); }
    VECTOR(*samples)[i] = ss;
    SPLICING_CHECK(splicing_matrix_init(ss, noiso, noSamples));
  }  
  
  SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));

  /* Initialize Psi(0) randomly */

  SPLICING_CHECK(splicing_drift_proposal_init((int) noiso, noChains,
					      psi, alpha, &sigma,
					      start, start_psi,
					      gff, (int) gene, readLength,
					      overHang, reads,
					      /*paired=*/ 0, 0, 0, 0, 0, 0));

  SPLICING_CHECK(splicing_drift_proposal_propose((int) noiso, noChains,
						 alpha, sigma, psi, alpha));
  
  /* Initialize assignments of reads */
  if (algorithm == SPLICING_ALGO_REASSIGN) {
    SPLICING_CHECK(splicing_reassign_samples(mymatch_matrix, &match_order,
					     psi, (int) noiso, noChains, &vass));
  }

  while (1) {

    for (m=0, rundata->noAccepted=0, rundata->noRejected=0; 
	 m < noIterations; 
	 m++) {
      
      SPLICING_CHECK(splicing_drift_proposal_propose((int) noiso, noChains,
						     alpha, sigma, 
						     psiNew, alphaNew));

      for (i = 0; i < noReplicates; i++) {

	splicing_matrix_int_t *vass1 = have_vass ? VECTOR(vass)[i] : 0;
	splicing_matrix_t *psiNew1 = VECTOR(*psiNew)[i];
	splicing_matrix_t *alphaNew1 = VECTOR(*alphaNew)[i];
	splicing_matrix_t *psi1 = VECTOR(*psi)[i];
	splicing_matrix_t *alpha1 = VECTOR(*alpha)[i];
	splicing_matrix_t *match_matrix = VECTOR(*mymatch_matrix)[i];
	splicing_vector_t *matches1 =
	  algorithm == SPLICING_ALGO_CLASSES ? VECTOR(matches)[i] : 0;

	SPLICING_CHECK(splicing_metropolis_hastings_ratio(algorithm,
	  vass1, VECTOR(noReads)[i], noChains, psiNew1, alphaNew1, psi1, alpha1,
	  sigma, (int) noiso, match_matrix, &assignmentMatrix, matches1,
	  &effisolen, hyperprior, &isoscores, m > 0 ? 1 : 0, &acceptP, &cJS, 
	  &pJS));
	

	for (j=0; j<noChains; j++) {
	  if (VECTOR(acceptP)[j] >= 1 || RNG_UNIF01() < VECTOR(acceptP)[j]) {
	    memcpy(&MATRIX(*psi1, 0, j), &MATRIX(*psiNew1, 0, j), 
		   noiso * sizeof(double));
	    memcpy(&MATRIX(*alpha1, 0, j), &MATRIX(*alphaNew1, 0, j),
		   (noiso - 1) * sizeof(double));
	    VECTOR(cJS)[j] = VECTOR(pJS)[j];
	    rundata->noAccepted ++;
	  } else {
	    rundata->noRejected ++;
	  }
	}
	
      }	/* i < noReplicates */
      
      if (m >= noBurnIn) {
	if (lagCounter == noLag - 1) {
	  for (i = 0; i < noReplicates; i++) {
	    splicing_matrix_t *samples1 = VECTOR(*samples)[i];
	    splicing_matrix_t *psi1 = VECTOR(*psi)[i];
	    memcpy(&MATRIX(*samples1, 0, noS), &MATRIX(*psi1, 0, 0), 
		   noChains * noiso * sizeof(double));
	    memcpy(VECTOR(*logLik)+noS, VECTOR(cJS), 
		   noChains * sizeof(double));
	    noS += noChains;
	    lagCounter = 0;
	  }
	} else {
	  lagCounter ++;
	}
      }
      
      if (algorithm == SPLICING_ALGO_REASSIGN) {
	SPLICING_CHECK(splicing_reassign_samples(mymatch_matrix, &match_order,
						 psi, (int) noiso, noChains, &vass));
      }
      
    } /* for m < noIterations */
    
    /* Should we stop? */
    switch (stop) {
    case SPLICING_MISO_STOP_FIXEDNO: 
      shouldstop = 1;
      break;
    case SPLICING_MISO_STOP_CONVERGENT_MEAN:
      if (maxIterations <= noIterations) { 
	shouldstop = 1; 
      } else {
	for (i = 0; i < noReplicates; i++) {
	  int shouldstop1 = 1;
	  splicing_matrix_t *cm = VECTOR(chainMeans)[i];
	  splicing_matrix_t *cv = VECTOR(chainVars)[i];
	  splicing_matrix_t *sam = VECTOR(*samples)[i];
	  SPLICING_CHECK(splicing_i_check_convergent_mean(cm, cv, sam,
							  &shouldstop1));
	  shouldstop = shouldstop && shouldstop1;
	}
      }
      break;
    }
    
    if (shouldstop) { break; }
    
    noS=0;
    noIterations = 3*noIterations - 2*noBurnIn;
    noBurnIn = m;
    noSamples = noChains * (noIterations - noBurnIn) / noLag;
    lagCounter = 0;

    for (i = 0; i < noReplicates; i++) {
      splicing_matrix_t *sam = VECTOR(*samples)[i];
      SPLICING_CHECK(splicing_matrix_resize(sam, noiso, noSamples));
    }
    SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));
  }
  
  splicing_vector_ptr_destroy(&chainVars);
  splicing_vector_ptr_destroy(&chainMeans);
  SPLICING_FINALLY_CLEAN(2);

  if (assignment) {
    splicing_vector_ptr_resize(assignment, noReplicates);
    splicing_vector_ptr_set_item_destructor(assignment,
      (splicing_finally_func_t *) splicing_vector_int_destroy_free);
    SPLICING_FINALLY(splicing_vector_ptr_destroy, assignment);
    /* This might not have been calculated, so we calculate it now,
       we could speed it up to only do it for a single chain, but
       it is done only once, so it does not matter much. */
    if (algorithm != SPLICING_ALGO_REASSIGN) {
      SPLICING_CHECK(splicing_reassign_samples(mymatch_matrix, &match_order,
					       psi, (int) noiso, noChains, &vass));
    }
    for (i = 0; i < noReplicates; i++) {
      int j, nor = VECTOR(noReads)[i];
      splicing_matrix_int_t *vass1 = VECTOR(vass)[i];
      splicing_vector_int_t *ass = splicing_Calloc(1, splicing_vector_int_t);
      if (!ass) { SPLICING_ERROR("No memory for assignment", SPLICING_ENOMEM); }
      VECTOR(*assignment)[i] = ass;
      SPLICING_CHECK(splicing_vector_int_init(ass, nor));
      for (j=0; j<nor; j++) {
	VECTOR(*ass)[j] = MATRIX(*vass1, j, 0);
      }
    }
  }

  if (algorithm == SPLICING_ALGO_CLASSES) {
    splicing_vector_ptr_destroy(&matches);
    splicing_matrix_destroy(&assignmentMatrix);
    SPLICING_FINALLY_CLEAN(2);
  }

  splicing_vector_destroy(&isoscores);
  splicing_vector_int_destroy(&effisolen);
  splicing_vector_ptr_destroy(&match_order);
  SPLICING_FINALLY_CLEAN(3);
  if (!match_matrix) {
    splicing_vector_ptr_destroy(mymatch_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  splicing_vector_ptr_destroy(&valphaNew);
  splicing_vector_ptr_destroy(&valpha);
  splicing_vector_ptr_destroy(&vpsiNew);
  splicing_vector_ptr_destroy(&vpsi);
  SPLICING_FINALLY_CLEAN(4);
  if (algorithm == SPLICING_ALGO_REASSIGN || assignment) {
    splicing_vector_ptr_destroy(&vass);
    SPLICING_FINALLY_CLEAN(1);
  }
  splicing_vector_destroy(&cJS);
  splicing_vector_destroy(&pJS);
  splicing_vector_destroy(&acceptP);
  SPLICING_FINALLY_CLEAN(3);

  /* We always return the same number of samples. */

  if (rundata->noSamples != noSamples) {
    splicing_vector_remove_section(logLik, 0, noSamples-rundata->noSamples);
    for (i = 0; i < noReplicates; i++) {
      splicing_matrix_t *sam = VECTOR(*samples)[i];
      splicing_matrix_remove_cols_section(sam, 0,
					  noSamples-rundata->noSamples);
    }
  }
  
  return 0;
}

int splicing_order_matches(const splicing_matrix_t *matches,
			   splicing_vector_int_t *order) {

  SPLICING_CHECK(splicing_matrix_order_cols(matches, order));
  return 0;
}

