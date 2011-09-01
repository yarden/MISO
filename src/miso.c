
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
	sumpsi += VECTOR(*psi)[j];					\
	VECTOR(validIso)[noValid] = j;					\
	VECTOR(cumsum)[noValid] = sumpsi;				\
	noValid++;							\
      }									\
    }									\
  } while (0)


int splicing_reassign_samples(const splicing_matrix_t *matches, 
			      const splicing_vector_int_t *match_order,
			      const splicing_vector_t *psi, 
			      int noiso, splicing_vector_int_t *result) {

  int noreads = splicing_matrix_ncol(matches);
  int i, w;
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

  SPLICING_CHECK(splicing_vector_int_resize(result, noreads));

  if (noreads == 0) { return 0; }  

  prev = curr = &MATRIX(*matches, 0, order[0]);
  CUMSUM();

  for (i=0; i<noreads; i++) {
    curr = &MATRIX(*matches, 0, order[i]);

    /* Maybe we need to update the cumulative sum */
    if (memcmp(prev, curr, sizeof(double)*noiso) != 0) { CUMSUM(); }
    
    if (noValid == 1) {
      VECTOR(*result)[order[i]] = VECTOR(validIso)[0];
    } else if (noValid == 2) { 
      rand = RNG_UNIF01() * sumpsi;
      w = (rand < VECTOR(cumsum)[0]) ? VECTOR(validIso)[0] : 
	VECTOR(validIso)[1];
      VECTOR(*result)[order[i]] = w;
    } else {
      /* Draw */
      rand = RNG_UNIF01() * sumpsi;
      /* TODO: Binary search for interval, if many classes */
      for (w=0; rand > VECTOR(cumsum)[w]; w++) ;
      VECTOR(*result)[order[i]] = VECTOR(validIso)[w];
    }

    prev=curr;
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
    score += VECTOR(logpsi)[ VECTOR(*assignment)[i] ];
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

int splicing_mvrnorm(const splicing_vector_t *mu, double sigma, 
		     splicing_vector_t *resalpha, int len) {
  int i;
  double sqrtsigma = len == 1 ? sigma : sqrt(sigma);

  SPLICING_CHECK(splicing_vector_resize(resalpha, len));

  for (i=0; i<len; i++) {
    VECTOR(*resalpha)[i] = VECTOR(*mu)[i] + sqrtsigma * RNG_NORMAL(0,1);
  }

  return 0;
}

int splicing_logit_inv(const splicing_vector_t *x, 
		       splicing_vector_t *res, int len) {
  int i;
  double sumexp=0.0;

  SPLICING_CHECK(splicing_vector_resize(res, len));

  for (i=0; i<len; i++) {
    sumexp += exp(VECTOR(*x)[i]);
  }
  sumexp += 1.0;
  
  for (i=0; i<len; i++) {
    VECTOR(*res)[i] = exp(VECTOR(*x)[i]) / sumexp;
  }
  
  return 0;
}

int splicing_score_joint(const splicing_vector_int_t *assignment,
			 int no_reads, const splicing_vector_t *psi, 
			 const splicing_vector_t *hyper, 
			 const splicing_vector_int_t *effisolen,
			 const splicing_vector_t *isoscores, 
			 double *score) {

  int i, noiso = splicing_vector_int_size(effisolen);
  double readProb = 0.0, assProb, psiProb;
  
  /* Scores the reads */
  for (i=0; i<no_reads; i++) {
    readProb += VECTOR(*isoscores)[ VECTOR(*assignment)[i] ];
  }
  
  /* Score isoforms */
  SPLICING_CHECK(splicing_score_iso(psi, noiso, assignment, no_reads, 
				    effisolen, &assProb));
  SPLICING_CHECK(splicing_ldirichlet(psi, hyper, noiso, &psiProb));

  *score = readProb + assProb + psiProb;
  return 0;
}

int splicing_drift_proposal(int mode, 
			    const splicing_vector_t *psi, 
			    const splicing_vector_t *alpha, 
			    double sigma, 
			    const splicing_vector_t *otherpsi, 
			    const splicing_vector_t *otheralpha, int noiso,
			    splicing_vector_t *respsi, 
			    splicing_vector_t *resalpha,
			    double *ressigma, double *resscore) {

  switch (mode) {
  case 0: 			/* init */
    {
      SPLICING_CHECK(splicing_vector_resize(respsi, noiso));
      SPLICING_CHECK(splicing_vector_resize(resalpha, noiso-1));
      if (noiso != 2) {
	int i;
	for (i=0; i<noiso; i++) {	
	  VECTOR(*respsi)[i] = 1.0/noiso; 
	}
	for (i=0; i<noiso-1; i++) { 
	  VECTOR(*resalpha)[i] = 1.0/(noiso-1);
	}
	*ressigma = 0.05;
      } else {
	VECTOR(*respsi)[0] = RNG_UNIF01();
	VECTOR(*respsi)[1] = 1 - VECTOR(*respsi)[0];
	VECTOR(*resalpha)[0] = 0.0;
	VECTOR(*resalpha)[1] = 0.0;
	*ressigma = 0.05;
      }
    }
    break;
  case 1:			/* propose */
    {
      int len=noiso-1;
      double sumpsi=0.0;
  
      SPLICING_CHECK(splicing_vector_reserve(respsi, len+1));
      SPLICING_CHECK(splicing_mvrnorm(alpha, sigma, resalpha, len));
      SPLICING_CHECK(splicing_logit_inv(resalpha, respsi, len));
      sumpsi = splicing_vector_sum(respsi);
      SPLICING_CHECK(splicing_vector_resize(respsi, len+1));
      VECTOR(*respsi)[len] = 1-sumpsi;
    }
    break;
  case 2: 			/* score */
    SPLICING_CHECK(splicing_mvplogisnorm(psi, otheralpha, sigma, noiso-1, 
					 resscore));
    break;
  }
  
  return 0;
}

int splicing_metropolis_hastings_ratio(const splicing_vector_int_t *ass,
				       int no_reads,
				       const splicing_vector_t *psiNew,
				       const splicing_vector_t *alphaNew,
				       const splicing_vector_t *psi, 
				       const splicing_vector_t *alpha,
				       double sigma,
				       int noiso, 
				       const splicing_vector_int_t *effisolen,
				       const splicing_vector_t *hyperp, 
				       const splicing_vector_t *isoscores,
				       int full, double *acceptP, 
				       double *pcJS, double *ppJS) {
  double pJS, cJS, ptoCS, ctoPS;

  SPLICING_CHECK(splicing_score_joint(ass, no_reads, psiNew, hyperp,
				      effisolen, isoscores, &pJS));
  SPLICING_CHECK(splicing_score_joint(ass, no_reads, psi, hyperp,
				      effisolen, isoscores, &cJS));
  
  SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 2, psi, alpha, sigma, 
					 psiNew, alphaNew, noiso, 0, 0, 0,
					 &ptoCS));
  SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 2, psiNew, alphaNew,
					 sigma, psi, alpha, noiso, 0, 0, 0,
					 &ctoPS));
  
  if (full) {
    *acceptP = exp(pJS + ptoCS - (cJS + ctoPS));
  } else {
    *acceptP = exp(pJS - cJS);
  }

  *pcJS = cJS;
  *ppJS = pJS;

  return 0;
}

int splicing_i_miso_classes(const splicing_matrix_t *match_matrix,
			    const splicing_vector_int_t *match_order,
			    splicing_matrix_t *class_templates,
			    splicing_vector_t *class_counts) {
  
  int noiso=splicing_matrix_nrow(match_matrix);

  if (splicing_matrix_size(match_matrix) == 0) { 

    /* Special case: no reads */
    splicing_matrix_resize(class_templates, noiso, 0);
    splicing_vector_resize(class_counts, 0);

  } else { 

    int i, noreads=splicing_vector_int_size(match_order);
    int lastclass=0;
    double *prev, *curr;
    int *order=VECTOR(*match_order);

    splicing_matrix_resize(class_templates, noiso, 1);
    splicing_vector_resize(class_counts, 1);
    
    memcpy(&MATRIX(*class_templates, 0, lastclass), 
	   &MATRIX(*match_matrix, 0, order[0]), 
	   sizeof(double) * noiso);
    prev = &MATRIX(*class_templates, 0, lastclass);
	   
    for (i=0; i<noreads; i++) {
      curr = &MATRIX(*match_matrix, 0, order[i]);
      if (memcmp(prev, curr, sizeof(double)*noiso) != 0) {
	SPLICING_CHECK(splicing_matrix_add_cols(class_templates, 1));
	SPLICING_CHECK(splicing_vector_push_back(class_counts, 0));
	lastclass++;
	prev = &MATRIX(*class_templates, 0, lastclass);
	memcpy(prev, curr, sizeof(double) * noiso);
      }
      VECTOR(*class_counts)[lastclass] += 1;
    }

  } 

  return 0;
}

int splicing_miso(const splicing_gff_t *gff, size_t gene,
		  const splicing_vector_int_t *position,
		  const char **cigarstr, int readLength, 
		  int noIterations, int noBurnIn, int noLag,
		  const splicing_vector_t *hyperp, 
		  splicing_matrix_t *samples, splicing_vector_t *logLik,
		  splicing_matrix_t *match_matrix, 
		  splicing_matrix_t *class_templates,
		  splicing_vector_t *class_counts,
		  splicing_miso_rundata_t *rundata) {

  double acceptP, cJS, pJS, sigma;
  int noReads = splicing_vector_int_size(position);
  splicing_vector_int_t ass;
  size_t noiso;
  splicing_vector_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;
  int noSamples = (noIterations - noBurnIn + 1) / noLag;
  int i, m, lagCounter=0, noS=0;
  splicing_matrix_t *mymatch_matrix=match_matrix, vmatch_matrix;
  splicing_vector_int_t match_order;
  splicing_vector_int_t effisolen;
  splicing_vector_t isoscores;

  if ( (class_templates ? 1 : 0) + (class_counts ? 1 : 0) == 1) {
    SPLICING_ERROR("Only one of `class_templates' and `class_counts' is "
		   "given", SPLICING_EINVAL);
  }

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));

  rundata->noIso=noiso;
  rundata->noIters=noIterations;
  rundata->noBurnIn=noBurnIn;
  rundata->noLag=noLag;
  rundata->noAccepted = rundata->noRejected = 0;

  SPLICING_CHECK(splicing_vector_int_init(&ass, noReads));
  SPLICING_FINALLY(splicing_vector_int_destroy, &ass);
  SPLICING_CHECK(splicing_vector_init(&vpsi, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &vpsi);
  SPLICING_CHECK(splicing_vector_init(&vpsiNew, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &vpsiNew);
  SPLICING_CHECK(splicing_vector_init(&valpha, noiso-1));
  SPLICING_FINALLY(splicing_vector_destroy, &valpha);
  SPLICING_CHECK(splicing_vector_init(&valphaNew, noiso-1));
  SPLICING_FINALLY(splicing_vector_destroy, &valphaNew);
  
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
				   mymatch_matrix));
  SPLICING_CHECK(splicing_order_matches(mymatch_matrix, &match_order));

  if (class_templates && class_counts) { 
    SPLICING_CHECK(splicing_i_miso_classes(mymatch_matrix, &match_order, 
					   class_templates, class_counts));
  }

  SPLICING_CHECK(splicing_vector_int_init(&effisolen, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &effisolen);
  SPLICING_CHECK(splicing_vector_init(&isoscores, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &isoscores);
  SPLICING_CHECK(splicing_gff_isolength_one(gff, gene, &effisolen));
  for (i=0; i<noiso; i++) { 
    int l=VECTOR(effisolen)[i]-readLength+1;
    VECTOR(effisolen)[i] = l > 0 ? l : 0;
    VECTOR(isoscores)[i] = -log((double) l);
  }

  SPLICING_CHECK(splicing_matrix_resize(samples, noiso, noSamples));
  SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));

  /* Initialize Psi(0) randomly */

  SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 0, 0, 0, 0, 0, 0, 
					 noiso, psi, alpha, &sigma, 0));
  SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 1, psi, alpha, sigma,
					 0, 0, noiso, psi, alpha, 0, 0));
  
  /* Initialize assignments of reads */  
  
  SPLICING_CHECK(splicing_reassign_samples(mymatch_matrix, &match_order, psi, 
					   noiso, &ass));
  
  /* foreach Iteration m=1, ..., M do */

  for (m=0; m < noIterations; m++) {

    SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 1, psi, alpha, sigma,
					   0, 0, noiso, psiNew, alphaNew, 0,
					   0));

    SPLICING_CHECK(splicing_metropolis_hastings_ratio(&ass, noReads, psiNew,
						      alphaNew, psi, alpha,
						      sigma, noiso, 
						      &effisolen, hyperp,
						      &isoscores, 
						      m > 0 ? 1 : 0, 
						      &acceptP, &cJS, &pJS));
    
    if (acceptP >= 1 || RNG_UNIF01() < acceptP) {
      splicing_vector_t *tmp;
      tmp=psi; psi=psiNew; psiNew=tmp;
      tmp=alpha; alpha=alphaNew; alphaNew=tmp;
      cJS = pJS;
      rundata->noAccepted ++;
    } else {
      rundata->noRejected ++;
    }
    
    if (m >= noBurnIn) {
      if (lagCounter == noLag - 1) {
	memcpy(&MATRIX(*samples, 0, noS), VECTOR(*psi), 
	       noiso * sizeof(double));
	VECTOR(*logLik)[noS] = cJS;
	noS++;
	lagCounter = 0;
      } else {
	lagCounter ++;
      }
    }
    
    SPLICING_CHECK(splicing_reassign_samples(mymatch_matrix, &match_order, 
					     psi, noiso, &ass));

  } /* for m < noIterations */

  splicing_vector_destroy(&isoscores);
  splicing_vector_int_destroy(&effisolen);
  splicing_vector_int_destroy(&match_order);
  SPLICING_FINALLY_CLEAN(3);
  if (!match_matrix) {
    splicing_matrix_destroy(mymatch_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  splicing_vector_destroy(&valphaNew);
  splicing_vector_destroy(&valpha);
  splicing_vector_destroy(&vpsiNew);
  splicing_vector_destroy(&vpsi);
  splicing_vector_int_destroy(&ass);
  SPLICING_FINALLY_CLEAN(5);

  return 0;
}

int splicing_order_matches(const splicing_matrix_t *matches,
			   splicing_vector_int_t *order) {

  SPLICING_CHECK(splicing_matrix_order_cols(matches, order));
  return 0;
}

