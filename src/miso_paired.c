
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
	sumpsi += VECTOR(*psi)[j] * (*cptr);				\
	VECTOR(validIso)[noValid] = j;					\
	VECTOR(cumsum)[noValid] = sumpsi;				\
	noValid++;							\
      }									\
    }									\
  } while (0)

int splicing_reassign_samples_paired(
			     const splicing_matrix_t *matches, 
			     const splicing_vector_int_t *match_order,
			     const splicing_vector_t *psi, 
			     int noiso, int insertStart, 
			     splicing_vector_int_t *result) {

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

int splicing_score_iso_paired(const splicing_vector_t *psi, int noiso, 
			      const splicing_vector_int_t *assignment, 
			      const splicing_vector_int_t *pisolen, 
			      const splicing_vector_t *assscores,
			      double *res) {

  int noreads=splicing_vector_int_size(assignment);
  int *isolen = VECTOR(*pisolen);
  double sum, maxpsieff, score;
  splicing_vector_t logpsi;
  int i;

  SPLICING_CHECK(splicing_vector_init(&logpsi, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &logpsi);

  /* Calculate the normalization factor */
  VECTOR(logpsi)[0] = log(VECTOR(*psi)[0]) + VECTOR(*assscores)[0];
  for (maxpsieff=VECTOR(logpsi)[0], i=1; i<noiso; i++) {
    VECTOR(logpsi)[i] = log(VECTOR(*psi)[i]) + VECTOR(*assscores)[i];
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

int splicing_score_joint_paired(const splicing_vector_int_t *assignment,
				const splicing_vector_t *psi, 
				const splicing_vector_t *hyper, 
				const splicing_vector_int_t *isolen,
				const splicing_matrix_t *isoscores, 
				const splicing_vector_t *assscores,
				const splicing_matrix_int_t *insertLength,
				int insertStart, double *score) {

  int no_reads=splicing_vector_int_size(assignment);
  int i, noiso = splicing_vector_int_size(isolen);
  double readProb = 0.0, assProb, psiProb;
  
  /* Scores the reads */
  for (i=0; i<no_reads; i++) {
    int ass=VECTOR(*assignment)[i];
    int inslen=MATRIX(*insertLength, ass, i);
    readProb += MATRIX(*isoscores, inslen-insertStart, ass);
  }
  
  /* Score isoforms */
  SPLICING_CHECK(splicing_score_iso_paired(psi, noiso, assignment,
					   isolen, assscores, &assProb));
  SPLICING_CHECK(splicing_ldirichlet(psi, hyper, noiso, &psiProb));

  *score = readProb + assProb + psiProb;
  return 0;
}

int splicing_metropolis_hastings_ratio_paired(
			      const splicing_vector_int_t *ass,
			      const splicing_vector_t *psiNew,
			      const splicing_vector_t *alphaNew,
			      const splicing_vector_t *psi, 
			      const splicing_vector_t *alpha,
			      double sigma, int noiso, 
			      const splicing_vector_int_t *isolen,
			      const splicing_vector_t *hyperp, 
			      const splicing_matrix_t *isoscores,
			      const splicing_vector_t *assscores,
			      const splicing_matrix_int_t *insertLength,
			      int insertStart, int full, double *acceptP, 
			      double *pcJS, double *ppJS) {
  
  double pJS, cJS, ptoCS, ctoPS;

  SPLICING_CHECK(splicing_score_joint_paired(ass, psiNew, hyperp,
					     isolen, isoscores, assscores, 
					     insertLength, insertStart,
					     &pJS));
  SPLICING_CHECK(splicing_score_joint_paired(ass, psi, hyperp,
					     isolen, isoscores, assscores,
					     insertLength, insertStart,
					     &cJS));
  
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

int splicing_miso_paired(const splicing_gff_t *gff, size_t gene,
			 const splicing_vector_int_t *position,
			 const char **cigarstr, int readLength, 
			 int noIterations, int noBurnIn, int noLag,
			 const splicing_vector_t *hyperp,
			 const splicing_vector_t *insertProb,
			 int insertStart, double normalMean, 
			 double normalVar, double numDevs,
			 splicing_matrix_t *samples, 
			 splicing_vector_t *logLik,
			 splicing_miso_rundata_t *rundata) {

  double acceptP, cJS, pJS, sigma;
  int noReads = splicing_vector_int_size(position)/2;
  splicing_vector_int_t ass;
  size_t noiso;
  splicing_vector_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;
  int noSamples = (noIterations - noBurnIn + 1) / noLag;
  int i, j, m, lagCounter=0, noS=0;
  splicing_matrix_t matches;
  splicing_vector_int_t match_order;
  splicing_vector_int_t isolen;
  splicing_matrix_t isoscores;
  splicing_matrix_int_t insertLength;
  splicing_vector_t assscores;
  int il;
  splicing_vector_t *myinsertProb=(splicing_vector_t*) insertProb,
    vinsertProb;

  if (!insertProb) { 
    myinsertProb=&vinsertProb;
    SPLICING_CHECK(splicing_vector_init(&vinsertProb, 0));
    SPLICING_FINALLY(splicing_vector_destroy, &vinsertProb);
    SPLICING_CHECK(splicing_normal_insert(normalMean, normalVar, numDevs,
					  myinsertProb, &insertStart));
    splicing_vector_scale(myinsertProb, 
			  1.0/splicing_vector_sum(myinsertProb));
  }

  il=splicing_vector_size(myinsertProb);

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

  SPLICING_CHECK(splicing_matrix_init(&matches, noiso, noReads));
  SPLICING_FINALLY(splicing_matrix_destroy, &matches);
  SPLICING_CHECK(splicing_vector_int_init(&match_order, noReads));
  SPLICING_FINALLY(splicing_vector_int_destroy, &match_order);
  SPLICING_CHECK(splicing_matrix_int_init(&insertLength, noiso, noReads));
  SPLICING_FINALLY(splicing_matrix_int_destroy, &insertLength);
  SPLICING_CHECK(splicing_matchIso_paired(gff, gene, position, cigarstr, 
					  readLength, myinsertProb, 
					  insertStart, normalMean, normalVar,
					  numDevs, &matches, &insertLength));
  SPLICING_CHECK(splicing_order_matches(&matches, &match_order));

  SPLICING_CHECK(splicing_vector_int_init(&isolen, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isolen);
  SPLICING_CHECK(splicing_gff_isolength_one(gff, gene, &isolen));
  SPLICING_CHECK(splicing_matrix_init(&isoscores, il, noiso));
  SPLICING_FINALLY(splicing_matrix_destroy, &isoscores);
  SPLICING_CHECK(splicing_vector_init(&assscores, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &assscores);
  for (j=0; j<il; j++) {
    double logprob=VECTOR(*myinsertProb)[j];
    for (i=0; i<noiso; i++) {
      double lp = VECTOR(isolen)[i] - insertStart - j - 2 * readLength + 1;
      MATRIX(isoscores, j, i) = -log(lp) + logprob;
      VECTOR(assscores)[i] += lp;
    }
  }
  for (i=0; i<noiso; i++) {
    VECTOR(assscores)[i] = log(VECTOR(assscores)[i]);
  }

  SPLICING_CHECK(splicing_matrix_resize(samples, noiso, noSamples));
  SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));

  /* Initialize Psi(0) randomly */

  SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 0, 0, 0, 0, 0, 0, noiso,
					 psi, alpha, &sigma, 0));
  SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 1, psi, alpha, sigma,
					 0, 0, noiso, psi, alpha, 0, 0));
  
  /* Initialize assignments of reads */  
  
  SPLICING_CHECK(splicing_reassign_samples_paired(&matches, &match_order, 
						  psi, noiso, insertStart,
						  &ass));
  
  /* foreach Iteration m=1, ..., M do */

  for (m=0; m < noIterations; m++) {

    SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 1, psi, alpha, sigma,
					   0, 0, noiso, psiNew, alphaNew,
					   0, 0));

    SPLICING_CHECK(splicing_metropolis_hastings_ratio_paired(&ass,
					     psiNew, alphaNew, psi, alpha,
					     sigma, noiso, &isolen, hyperp, 
					     &isoscores, &assscores,
					     &insertLength, insertStart,
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
    
    SPLICING_CHECK(splicing_reassign_samples_paired(&matches, &match_order,
						    psi, noiso, insertStart,
						    &ass));

  } /* for m < noIterations */

  splicing_vector_destroy(&assscores);
  splicing_matrix_destroy(&isoscores);
  splicing_vector_int_destroy(&isolen);
  splicing_matrix_int_destroy(&insertLength);
  splicing_vector_int_destroy(&match_order);
  splicing_matrix_destroy(&matches);
  splicing_vector_destroy(&valphaNew);
  splicing_vector_destroy(&valpha);
  splicing_vector_destroy(&vpsiNew);
  splicing_vector_destroy(&vpsi);
  splicing_vector_int_destroy(&ass);
  SPLICING_FINALLY_CLEAN(11);

  if (!insertProb) { 
    splicing_vector_destroy(&vinsertProb);
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}
