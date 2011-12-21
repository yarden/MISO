
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "splicing_matrix.h"
#include "splicing_error.h"
#include "splicing.h"
#include "splicing_random.h"

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


/* TODO: we could actually speed this up, by storing the 
   indices where the cumulative sum vector must be updated. 
   These indices do not change at all. With this we can spare
   the memcmp operations. */

int splicing_reassign_samples(const splicing_matrix_t *matches, 
			      const splicing_vector_int_t *match_order,
			      const splicing_matrix_t *psi, 
			      int noiso, int noChains, 
			      splicing_matrix_int_t *result) {

  int noreads = splicing_matrix_ncol(matches);
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

int splicing_mvrnorm(const splicing_matrix_t *mu, double sigma, 
		     splicing_matrix_t *resalpha, int len) {
  int i, j;
  int noChains=splicing_matrix_ncol(mu);
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

int splicing_score_joint(const splicing_matrix_int_t *assignment,
			 int no_reads, int noChains, 
			 const splicing_matrix_t *psi, 
			 const splicing_vector_t *hyper, 
			 const splicing_vector_int_t *effisolen,
			 const splicing_vector_t *isoscores, 
			 splicing_vector_t *score) {

  int i, j, noiso = splicing_vector_int_size(effisolen);

  SPLICING_CHECK(splicing_vector_resize(score, noChains));

  for (j=0; j<noChains; j++) {
    double readProb = 0.0, assProb, psiProb;
    splicing_vector_t tmp;
    splicing_vector_int_t tmp2;
    splicing_vector_view(&tmp, &MATRIX(*psi, 0, j), noiso);
    splicing_vector_int_view(&tmp2, &MATRIX(*assignment, 0, j), no_reads);

    /* Scores the reads */
    for (i=0; i<no_reads; i++) {
      if (MATRIX(*assignment, i, j) != -1) {
	readProb += VECTOR(*isoscores)[ MATRIX(*assignment, i, j) ];
      }
    }
    
    /* Score isoforms */
    SPLICING_CHECK(splicing_score_iso(&tmp, noiso, &tmp2, no_reads, 
				      effisolen, &assProb));
    SPLICING_CHECK(splicing_ldirichlet(&tmp, hyper, noiso, &psiProb));
    
    VECTOR(*score)[j] = readProb + assProb + psiProb;
  }

  return 0;
}

int splicing_rng_get_dirichlet(splicing_rng_t *rng, 
			       const splicing_vector_t *alpha, 
			       splicing_vector_t *result) {

  int i, l=splicing_vector_size(alpha);
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

#define SIGMA (0.2/noiso/noiso)

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
	MATRIX(*resalpha, 1, j) = 0.0;
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
				       /*match_matrix=*/ 0, 
				       /*assignment_matrix=*/ 0,
				       &tmp));
      } else {
	SPLICING_CHECK(splicing_solve_gene_paired(gff, gene, readLength,
				       overHang, position,
				       cigarstr, fragmentProb, fragmentStart,
				       normalMean, normalVar, numDevs,
				       /*match_matrix=*/ 0,
				       /*assignment_matrix=*/ 0, &tmp));
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

int splicing_drift_proposal_propose(int noiso, int noChains, 
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

  SPLICING_CHECK(splicing_score_joint(ass, no_reads, noChains, psiNew, 
				      hyperp, effisolen, isoscores, ppJS));
  SPLICING_CHECK(splicing_score_joint(ass, no_reads, noChains, psi, hyperp,
				      effisolen, isoscores, pcJS));
  
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
  int noiso = splicing_matrix_nrow(chainMeans);
  int noChains = splicing_matrix_ncol(chainMeans);
  int noSamples = splicing_matrix_ncol(samples);

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

int splicing_miso(const splicing_gff_t *gff, size_t gene,
		  const splicing_vector_int_t *position,
		  const char **cigarstr, int readLength, int overHang,
		  int noChains, int noIterations, 
		  int maxIterations, int noBurnIn, int noLag,
		  const splicing_vector_t *hyperp, 
		  splicing_miso_start_t start,
		  splicing_miso_stop_t stop,
		  const splicing_matrix_t *start_psi,
		  splicing_matrix_t *samples, splicing_vector_t *logLik,
		  splicing_matrix_t *match_matrix, 
		  splicing_matrix_t *class_templates,
		  splicing_vector_t *class_counts,
		  splicing_vector_int_t *assignment,
		  splicing_miso_rundata_t *rundata) {

  splicing_vector_t acceptP, cJS, pJS;
  double sigma;
  int noReads = splicing_vector_int_size(position);
  splicing_matrix_int_t vass;
  size_t noiso;
  splicing_matrix_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;
  int noSamples = noChains * (noIterations - noBurnIn) / noLag;  
  int i, j, m=0, lagCounter=0, noS=0;
  splicing_matrix_t *mymatch_matrix=match_matrix, vmatch_matrix;
  splicing_vector_int_t match_order;
  splicing_vector_int_t effisolen;
  splicing_vector_t isoscores;
  splicing_vector_int_t noexons;
  int shouldstop=0;
  splicing_matrix_t chainMeans, chainVars;

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

  if (splicing_vector_size(hyperp) != noiso) { 
    SPLICING_ERROR("Invalid hyperparameter vector length", 
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

  rundata->noIso=noiso;
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

  SPLICING_CHECK(splicing_matrix_int_init(&vass, noReads, noChains));
  SPLICING_FINALLY(splicing_matrix_int_destroy, &vass);
  SPLICING_CHECK(splicing_matrix_init(&vpsi, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &vpsi);
  SPLICING_CHECK(splicing_matrix_init(&vpsiNew, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &vpsiNew);
  SPLICING_CHECK(splicing_matrix_init(&valpha, noiso-1, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &valpha);
  SPLICING_CHECK(splicing_matrix_init(&valphaNew, noiso-1, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &valphaNew);
  
  if (match_matrix) { 
    SPLICING_CHECK(splicing_matrix_resize(match_matrix, noiso, noReads));
  } else {
    mymatch_matrix=&vmatch_matrix;
    SPLICING_CHECK(splicing_matrix_init(mymatch_matrix, noiso, noReads));
    SPLICING_FINALLY(splicing_matrix_destroy, mymatch_matrix);
  }
  SPLICING_CHECK(splicing_vector_int_init(&match_order, noReads));
  SPLICING_FINALLY(splicing_vector_int_destroy, &match_order);
  SPLICING_CHECK(splicing_matchIso(gff, gene, position, cigarstr, 
				   overHang, readLength, mymatch_matrix));
  SPLICING_CHECK(splicing_order_matches(mymatch_matrix, &match_order));

  if (class_templates && class_counts) { 
    SPLICING_CHECK(splicing_i_miso_classes(mymatch_matrix, &match_order, 
					   class_templates, class_counts, 
					   /*bin_class_templates=*/ 0,
					   /*bin_class_counts=*/ 0));
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

  SPLICING_CHECK(splicing_matrix_init(&chainMeans, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &chainMeans);
  SPLICING_CHECK(splicing_matrix_init(&chainVars, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &chainVars);
  
  SPLICING_CHECK(splicing_matrix_resize(samples, noiso, noSamples));
  SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));

  /* Initialize Psi(0) randomly */

  SPLICING_CHECK(splicing_drift_proposal_init(noiso, noChains, 
					      psi, alpha, &sigma,
					      start, start_psi,
					      gff, gene, readLength,
					      overHang, position, cigarstr,
					      /*paired=*/ 0, 0, 0, 0, 0, 0));

  SPLICING_CHECK(splicing_drift_proposal_propose(noiso, noChains, 
						 alpha, sigma, psi, alpha));
  
  /* Initialize assignments of reads */  
  
  SPLICING_CHECK(splicing_reassign_samples(mymatch_matrix, &match_order,
					   psi, noiso, noChains, &vass));

  while (1) {

    for (m=0, rundata->noAccepted=0, rundata->noRejected=0; 
	 m < noIterations; 
	 m++) {
      
      SPLICING_CHECK(splicing_drift_proposal_propose(noiso, noChains, 
						     alpha, sigma, 
						     psiNew, alphaNew));

      SPLICING_CHECK(splicing_metropolis_hastings_ratio(&vass, noReads, 
							noChains, psiNew,
							alphaNew, psi, alpha,
							sigma, noiso, 
							&effisolen, hyperp,
							&isoscores, 
							m > 0 ? 1 : 0, 
							&acceptP, &cJS, 
							&pJS));

      for (j=0; j<noChains; j++) {
	if (VECTOR(acceptP)[j] >= 1 || RNG_UNIF01() < VECTOR(acceptP)[j]) {
	  memcpy(&MATRIX(*psi, 0, j), &MATRIX(*psiNew, 0, j), 
		 noiso * sizeof(double));
	  memcpy(&MATRIX(*alpha, 0, j), &MATRIX(*alphaNew, 0, j),
		 (noiso - 1) * sizeof(double));
	  VECTOR(cJS)[j] = VECTOR(pJS)[j];
	  rundata->noAccepted ++;
	} else {
	  rundata->noRejected ++;
	}
      }
      
      if (m >= noBurnIn) {
	if (lagCounter == noLag - 1) {
	  memcpy(&MATRIX(*samples, 0, noS), &MATRIX(*psi, 0, 0), 
		 noChains * noiso * sizeof(double));
	  memcpy(VECTOR(*logLik)+noS, VECTOR(cJS), 
		 noChains * sizeof(double));
	  noS += noChains;
	  lagCounter = 0;
	} else {
	  lagCounter ++;
	}
      }
      
      SPLICING_CHECK(splicing_reassign_samples(mymatch_matrix, &match_order, 
					       psi, noiso, noChains, &vass));
      
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
	SPLICING_CHECK(splicing_i_check_convergent_mean(&chainMeans, 
							&chainVars, 
							samples, 
							&shouldstop));
      }
      break;
    }
    
    if (shouldstop) { break; }
    
    noS=0;
    noIterations = 3*noIterations - 2*noBurnIn;
    noBurnIn = m;
    rundata->noSamples = noSamples = 
      noChains * (noIterations - noBurnIn) / noLag;
    lagCounter = 0;

    SPLICING_CHECK(splicing_matrix_resize(samples, noiso, noSamples));
    SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));
  }
  
  splicing_matrix_destroy(&chainVars);
  splicing_matrix_destroy(&chainMeans);
  SPLICING_FINALLY_CLEAN(2);

  if (assignment) {
    SPLICING_CHECK(splicing_vector_int_resize(assignment, noReads));
    for (i=0; i<noReads; i++) {
      VECTOR(*assignment)[i] = MATRIX(vass, i, 0);
    }
  }

  splicing_vector_destroy(&isoscores);
  splicing_vector_int_destroy(&effisolen);
  splicing_vector_int_destroy(&match_order);
  SPLICING_FINALLY_CLEAN(3);
  if (!match_matrix) {
    splicing_matrix_destroy(mymatch_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  splicing_matrix_destroy(&valphaNew);
  splicing_matrix_destroy(&valpha);
  splicing_matrix_destroy(&vpsiNew);
  splicing_matrix_destroy(&vpsi);
  splicing_matrix_int_destroy(&vass);
  splicing_vector_destroy(&cJS);
  splicing_vector_destroy(&pJS);
  splicing_vector_destroy(&acceptP);
  SPLICING_FINALLY_CLEAN(8);
  
  return 0;
}

int splicing_order_matches(const splicing_matrix_t *matches,
			   splicing_vector_int_t *order) {

  SPLICING_CHECK(splicing_matrix_order_cols(matches, order));
  return 0;
}

