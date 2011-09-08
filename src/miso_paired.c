
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
			     int noiso, int fragmentStart, 
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
    
    if (noValid == 0) {
      VECTOR(*result)[order[i]] = -1;
    } else if (noValid == 1) {
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
    if (VECTOR(*assignment)[i] != -1) {
      score += VECTOR(logpsi)[ VECTOR(*assignment)[i] ];
    }
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
				const splicing_matrix_int_t *fragmentLength,
				int fragmentStart, double *score) {

  int no_reads=splicing_vector_int_size(assignment);
  int i, noiso = splicing_vector_int_size(isolen);
  double readProb = 0.0, assProb, psiProb;
  
  /* Scores the reads */
  for (i=0; i<no_reads; i++) {
    int ass=VECTOR(*assignment)[i];
    if (ass != -1) {
      int fraglen=MATRIX(*fragmentLength, ass, i);
      readProb += MATRIX(*isoscores, fraglen-fragmentStart, ass);
    }
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
			      const splicing_matrix_int_t *fragmentLength,
			      int fragmentStart, int full, double *acceptP, 
			      double *pcJS, double *ppJS) {
  
  double pJS, cJS, ptoCS, ctoPS;

  SPLICING_CHECK(splicing_score_joint_paired(ass, psiNew, hyperp,
					     isolen, isoscores, assscores, 
					     fragmentLength, fragmentStart,
					     &pJS));
  SPLICING_CHECK(splicing_score_joint_paired(ass, psi, hyperp,
					     isolen, isoscores, assscores,
					     fragmentLength, fragmentStart,
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
			 const splicing_vector_t *fragmentProb,
			 int fragmentStart, double normalMean, 
			 double normalVar, double numDevs,
			 splicing_matrix_t *samples, 
			 splicing_vector_t *logLik,
			 splicing_matrix_t *match_matrix, 
			 splicing_matrix_t *class_templates,
			 splicing_vector_t *class_counts,
			 splicing_vector_int_t *assignment,
			 splicing_miso_rundata_t *rundata) {

  double acceptP, cJS, pJS, sigma;
  int noReads = splicing_vector_int_size(position)/2;
  splicing_vector_int_t *myass=assignment, vass;
  size_t noiso;
  splicing_vector_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;
  int noSamples = (noIterations - noBurnIn + 1) / noLag;
  int i, j, m, lagCounter=0, noS=0;
  splicing_matrix_t *mymatch_matrix=match_matrix, vmatch_matrix;
  splicing_vector_int_t match_order;
  splicing_vector_int_t isolen;
  splicing_matrix_t isoscores;
  splicing_matrix_int_t fragmentLength;
  splicing_vector_t assscores;
  int il;
  splicing_vector_t *myfragmentProb=(splicing_vector_t*) fragmentProb,
    vfragmentProb;

  if ( (class_templates ? 1 : 0) + (class_counts ? 1 : 0) == 1) {
    SPLICING_ERROR("Only one of `class_templates' and `class_counts' is "
		   "given", SPLICING_EINVAL);
  }

  if (!fragmentProb) { 
    myfragmentProb=&vfragmentProb;
    SPLICING_CHECK(splicing_vector_init(&vfragmentProb, 0));
    SPLICING_FINALLY(splicing_vector_destroy, &vfragmentProb);
    SPLICING_CHECK(splicing_normal_fragment(normalMean, normalVar, numDevs,
					    2*readLength, myfragmentProb,
					    &fragmentStart));
    splicing_vector_scale(myfragmentProb, 
			  1.0/splicing_vector_sum(myfragmentProb));
  }

  il=splicing_vector_size(myfragmentProb);

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));

  rundata->noIso=noiso;
  rundata->noIters=noIterations;
  rundata->noBurnIn=noBurnIn;
  rundata->noLag=noLag;
  rundata->noAccepted = rundata->noRejected = 0;

  if (assignment) { 
    SPLICING_CHECK(splicing_vector_int_resize(myass, noReads));
    splicing_vector_int_null(myass);
  } else {
    myass=&vass;
    SPLICING_CHECK(splicing_vector_int_init(myass, noReads));
    SPLICING_FINALLY(splicing_vector_int_destroy, myass);
  }
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
  SPLICING_CHECK(splicing_matrix_int_init(&fragmentLength, noiso, noReads));
  SPLICING_FINALLY(splicing_matrix_int_destroy, &fragmentLength);
  SPLICING_CHECK(splicing_matchIso_paired(gff, gene, position, cigarstr, 
					  readLength, myfragmentProb, 
					  fragmentStart, normalMean, 
					  normalVar, numDevs, mymatch_matrix,
					  &fragmentLength));
  SPLICING_CHECK(splicing_order_matches(mymatch_matrix, &match_order));

  if (class_templates && class_counts) { 
    SPLICING_CHECK(splicing_i_miso_classes(mymatch_matrix, &match_order, 
					   class_templates, class_counts));
  }

  SPLICING_CHECK(splicing_vector_int_init(&isolen, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isolen);
  SPLICING_CHECK(splicing_gff_isolength_one(gff, gene, &isolen));
  SPLICING_CHECK(splicing_matrix_init(&isoscores, il, noiso));
  SPLICING_FINALLY(splicing_matrix_destroy, &isoscores);
  SPLICING_CHECK(splicing_vector_init(&assscores, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &assscores);
  for (j=0; j<il; j++) {
    double logprob=VECTOR(*myfragmentProb)[j];
    for (i=0; i<noiso; i++) {
      double lp = VECTOR(isolen)[i] - fragmentStart - j + 1;
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
  
  SPLICING_CHECK(splicing_reassign_samples_paired(mymatch_matrix,
						  &match_order, 
						  psi, noiso, fragmentStart,
						  myass));
  
  /* foreach Iteration m=1, ..., M do */

  for (m=0; m < noIterations; m++) {

    SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 1, psi, alpha, sigma,
					   0, 0, noiso, psiNew, alphaNew,
					   0, 0));

    SPLICING_CHECK(splicing_metropolis_hastings_ratio_paired(myass,
					     psiNew, alphaNew, psi, alpha,
					     sigma, noiso, &isolen, hyperp, 
					     &isoscores, &assscores,
					     &fragmentLength, fragmentStart,
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
    
    SPLICING_CHECK(splicing_reassign_samples_paired(mymatch_matrix,
						    &match_order,
						    psi, noiso, fragmentStart,
						    myass));

  } /* for m < noIterations */

  splicing_vector_destroy(&assscores);
  splicing_matrix_destroy(&isoscores);
  splicing_vector_int_destroy(&isolen);
  splicing_matrix_int_destroy(&fragmentLength);
  splicing_vector_int_destroy(&match_order);
  SPLICING_FINALLY_CLEAN(5);
  if (!match_matrix) {
    splicing_matrix_destroy(mymatch_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  splicing_vector_destroy(&valphaNew);
  splicing_vector_destroy(&valpha);
  splicing_vector_destroy(&vpsiNew);
  splicing_vector_destroy(&vpsi);
  SPLICING_FINALLY_CLEAN(4);
  if (!assignment) { 
    splicing_vector_int_destroy(myass);
    SPLICING_FINALLY_CLEAN(1);
  }

  if (!fragmentProb) { 
    splicing_vector_destroy(&vfragmentProb);
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}

int splicing_miso_paired_trinity(const splicing_matrix_t *match_matrix,
				 const splicing_matrix_int_t *fragmentLength,
				 const splicing_vector_int_t *isolen,
				 int readLength, int noIterations, 
				 int noBurnIn, int noLag, 
				 const splicing_vector_t *hyperp,
				 const splicing_vector_t *fragmentProb,
				 int fragmentStart, double normalMean,
				 double normalVar, double numDevs,
				 splicing_matrix_t *samples,
				 splicing_vector_t *logLik,
				 splicing_matrix_t *class_templates,
				 splicing_vector_t *class_counts,
				 splicing_vector_int_t *assignment,
				 splicing_miso_rundata_t *rundata) {
  
  double acceptP, cJS, pJS, sigma;
  int noiso = splicing_matrix_nrow(match_matrix);
  int noReads = splicing_matrix_ncol(match_matrix);
  splicing_vector_int_t *myass=assignment, vass;
  splicing_vector_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;
  int noSamples = (noIterations - noBurnIn + 1) / noLag;
  int i, j, m, lagCounter=0, noS=0;
  splicing_vector_int_t match_order;
  splicing_matrix_t isoscores;
  splicing_vector_t assscores;
  int il;
  splicing_vector_t *myfragmentProb=(splicing_vector_t*) fragmentProb,
    vfragmentProb;

  if ( (class_templates ? 1 : 0) + (class_counts ? 1 : 0) == 1) {
    SPLICING_ERROR("Only one of `class_templates' and `class_counts' is "
		   "given", SPLICING_EINVAL);
  }

  if (!fragmentProb) { 
    myfragmentProb=&vfragmentProb;
    SPLICING_CHECK(splicing_vector_init(&vfragmentProb, 0));
    SPLICING_FINALLY(splicing_vector_destroy, &vfragmentProb);
    SPLICING_CHECK(splicing_normal_fragment(normalMean, normalVar, numDevs,
					    2*readLength, myfragmentProb,
					    &fragmentStart));
    splicing_vector_scale(myfragmentProb, 
			  1.0/splicing_vector_sum(myfragmentProb));
  }

  il=splicing_vector_size(myfragmentProb);

  rundata->noIso=noiso;
  rundata->noIters=noIterations;
  rundata->noBurnIn=noBurnIn;
  rundata->noLag=noLag;
  rundata->noAccepted = rundata->noRejected = 0;

  if (assignment) { 
    SPLICING_CHECK(splicing_vector_int_resize(myass, noReads));
    splicing_vector_int_null(myass);
  } else {
    myass=&vass;
    SPLICING_CHECK(splicing_vector_int_init(myass, noReads));
    SPLICING_FINALLY(splicing_vector_int_destroy, myass);
  }
  SPLICING_CHECK(splicing_vector_init(&vpsi, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &vpsi);
  SPLICING_CHECK(splicing_vector_init(&vpsiNew, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &vpsiNew);
  SPLICING_CHECK(splicing_vector_init(&valpha, noiso-1));
  SPLICING_FINALLY(splicing_vector_destroy, &valpha);
  SPLICING_CHECK(splicing_vector_init(&valphaNew, noiso-1));
  SPLICING_FINALLY(splicing_vector_destroy, &valphaNew);

  SPLICING_CHECK(splicing_vector_int_init(&match_order, noReads));
  SPLICING_FINALLY(splicing_vector_int_destroy, &match_order);
  SPLICING_CHECK(splicing_order_matches(match_matrix, &match_order));
  
  if (class_templates && class_counts) { 
    SPLICING_CHECK(splicing_i_miso_classes(match_matrix, &match_order, 
					   class_templates, class_counts));
  }

  SPLICING_CHECK(splicing_matrix_init(&isoscores, il, noiso));
  SPLICING_FINALLY(splicing_matrix_destroy, &isoscores);
  SPLICING_CHECK(splicing_vector_init(&assscores, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &assscores);
  for (j=0; j<il; j++) {
    double logprob=VECTOR(*myfragmentProb)[j];
    for (i=0; i<noiso; i++) {
      double lp = VECTOR(*isolen)[i] - fragmentStart - j + 1;
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
  
  SPLICING_CHECK(splicing_reassign_samples_paired(match_matrix,
						  &match_order, 
						  psi, noiso, fragmentStart,
						  myass));

  /* foreach Iteration m=1, ..., M do */

  for (m=0; m < noIterations; m++) {

    SPLICING_CHECK(splicing_drift_proposal(/* mode= */ 1, psi, alpha, sigma,
					   0, 0, noiso, psiNew, alphaNew,
					   0, 0));

    SPLICING_CHECK(splicing_metropolis_hastings_ratio_paired(myass,
					     psiNew, alphaNew, psi, alpha,
					     sigma, noiso, isolen, hyperp, 
					     &isoscores, &assscores,
					     fragmentLength, fragmentStart,
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
    
    SPLICING_CHECK(splicing_reassign_samples_paired(match_matrix,
						    &match_order,
						    psi, noiso, fragmentStart,
						    myass));

  } /* for m < noIterations */

  splicing_vector_destroy(&assscores);
  splicing_matrix_destroy(&isoscores);
  splicing_vector_int_destroy(&match_order);
  splicing_vector_destroy(&valphaNew);
  splicing_vector_destroy(&valpha);
  splicing_vector_destroy(&vpsiNew);
  splicing_vector_destroy(&vpsi);
  SPLICING_FINALLY_CLEAN(7);
  if (!assignment) { 
    splicing_vector_int_destroy(myass);
    SPLICING_FINALLY_CLEAN(1);
  }

  if (!fragmentProb) { 
    splicing_vector_destroy(&vfragmentProb);
    SPLICING_FINALLY_CLEAN(1);
  }  
  
  return 0;
}
