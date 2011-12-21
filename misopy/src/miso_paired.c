
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
	sumpsi += MATRIX(*psi, j, k) * (*cptr);				\
	VECTOR(validIso)[noValid] = j;					\
	VECTOR(cumsum)[noValid] = sumpsi;				\
	noValid++;							\
      }									\
    }									\
  } while (0)

int splicing_reassign_samples_paired(
			     const splicing_matrix_t *matches, 
			     const splicing_vector_int_t *match_order,
			     const splicing_matrix_t *psi, 
			     int noiso, int noChains, int fragmentStart, 
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

int splicing_score_joint_paired(const splicing_matrix_int_t *assignment,
				int no_reads, int noChains,
				const splicing_matrix_t *psi, 
				const splicing_vector_t *hyper, 
				const splicing_vector_int_t *isolen,
				const splicing_matrix_t *isoscores, 
				const splicing_vector_t *assscores,
				const splicing_matrix_int_t *fragmentLength,
				int fragmentStart, 
				splicing_vector_t *score) {

  int i, j, noiso = splicing_vector_int_size(isolen);

  SPLICING_CHECK(splicing_vector_resize(score, noChains));

  /* Scores the reads */
  for (j=0; j<noChains; j++) {
    double readProb = 0.0, assProb, psiProb;
    splicing_vector_t tmp;
    splicing_vector_int_t tmp2;

    splicing_vector_view(&tmp, &MATRIX(*psi, 0, j), noiso);
    splicing_vector_int_view(&tmp2, &MATRIX(*assignment, 0, j), no_reads);

    for (i=0; i<no_reads; i++) {
      int ass=MATRIX(*assignment, i, j);
      if (ass != -1) {
	int fraglen=MATRIX(*fragmentLength, ass, i);
	readProb += MATRIX(*isoscores, fraglen-fragmentStart, ass);
      }
    }
    
    /* Score isoforms */
    SPLICING_CHECK(splicing_score_iso_paired(&tmp, noiso, &tmp2,
					     isolen, assscores, &assProb));
    SPLICING_CHECK(splicing_ldirichlet(&tmp, hyper, noiso, &psiProb));
    
    VECTOR(*score)[j] = readProb + assProb + psiProb;
  }
  
  return 0;
}

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

  SPLICING_CHECK(splicing_score_joint_paired(ass, no_reads, noChains, 
					     psiNew, hyperp,
					     isolen, isoscores, assscores, 
					     fragmentLength, fragmentStart,
					     ppJS));
  SPLICING_CHECK(splicing_score_joint_paired(ass, no_reads, noChains, 
					     psi, hyperp,
					     isolen, isoscores, assscores,
					     fragmentLength, fragmentStart,
					     pcJS));
  
  SPLICING_CHECK(splicing_drift_proposal_score(noiso, noChains, psi, 
					       alphaNew, sigma, &ptoCS));
  SPLICING_CHECK(splicing_drift_proposal_score(noiso, noChains, psiNew, 
					       alpha, sigma, &ctoPS));
  
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

int splicing_miso_paired(const splicing_gff_t *gff, size_t gene,
			 const splicing_vector_int_t *position,
			 const char **cigarstr, int readLength, int overHang,
			 int noChains, int noIterations, 
			 int maxIterations, int noBurnIn, int noLag,
			 const splicing_vector_t *hyperp, 
			 splicing_miso_start_t start, 
			 splicing_miso_stop_t stop,
			 const splicing_matrix_t *start_psi,
			 const splicing_vector_t *fragmentProb,
			 int fragmentStart, double normalMean, 
			 double normalVar, double numDevs,
			 splicing_matrix_t *samples, 
			 splicing_vector_t *logLik,
			 splicing_matrix_t *match_matrix, 
			 splicing_matrix_t *class_templates,
			 splicing_vector_t *class_counts,
			 splicing_matrix_t *bin_class_templates,
			 splicing_vector_t *bin_class_counts,
			 splicing_vector_int_t *assignment,
			 splicing_miso_rundata_t *rundata) {

  splicing_vector_t acceptP, cJS, pJS;
  double sigma;
  int noReads = splicing_vector_int_size(position)/2;
  splicing_matrix_int_t vass;
  size_t noiso;
  splicing_matrix_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;
  int noSamples = noChains * (noIterations - noBurnIn) / noLag;
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
  if ( (bin_class_templates ? 1 : 0) + (bin_class_counts ? 1 : 0) == 1) {
    SPLICING_ERROR("Only one of `bin_class_templates' and "
		   "`bin_class_counts' is given", SPLICING_EINVAL);
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

  if (splicing_vector_size(hyperp) != noiso) { 
    SPLICING_ERROR("Invalid hyperparameter vector length", 
		   SPLICING_EINVAL);
  }

  if (overHang==0) { overHang=1; }
  if (overHang < 1 || overHang >= readLength / 2) {
    SPLICING_ERROR("Overhang length invalid. Must be between 0 and "
		   "readLength/2", SPLICING_EINVAL);
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
  SPLICING_CHECK(splicing_matrix_int_init(&fragmentLength, noiso, noReads));
  SPLICING_FINALLY(splicing_matrix_int_destroy, &fragmentLength);
  SPLICING_CHECK(splicing_matchIso_paired(gff, gene, position, cigarstr, 
					  readLength, overHang, 
					  myfragmentProb, 
					  fragmentStart, normalMean, 
					  normalVar, numDevs, mymatch_matrix,
					  &fragmentLength));
  SPLICING_CHECK(splicing_order_matches(mymatch_matrix, &match_order));

  if (class_templates || bin_class_templates) {
    SPLICING_CHECK(splicing_i_miso_classes(mymatch_matrix, &match_order, 
					   class_templates, class_counts,
					   bin_class_templates, 
					   bin_class_counts));
  }

  SPLICING_CHECK(splicing_vector_int_init(&isolen, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isolen);
  SPLICING_CHECK(splicing_gff_isolength_one(gff, gene, &isolen));
  SPLICING_CHECK(splicing_matrix_init(&isoscores, il, noiso));
  SPLICING_FINALLY(splicing_matrix_destroy, &isoscores);
  SPLICING_CHECK(splicing_vector_init(&assscores, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &assscores);
  SPLICING_CHECK(splicing_vector_int_init(&noexons, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &noexons);
  SPLICING_CHECK(splicing_gff_noexons_one(gff, gene, &noexons));
  for (j=0; j<il; j++) {
    double logprob=VECTOR(*myfragmentProb)[j];
    for (i=0; i<noiso; i++) {
      int nox=VECTOR(noexons)[i];      
      /* The following is only approximate if there are some short exons
	 and overHang is not one */
      double lp = VECTOR(isolen)[i] - fragmentStart - j + 1 - 
	2 * (nox-1) * (overHang-1);
      MATRIX(isoscores, j, i) = -log(lp) + logprob;
      if (lp > 0) { VECTOR(assscores)[i] += lp; }
    }
  }
  splicing_vector_int_destroy(&noexons);
  SPLICING_FINALLY_CLEAN(1);
  for (i=0; i<noiso; i++) {
    VECTOR(assscores)[i] = log(VECTOR(assscores)[i]);
  }

  SPLICING_CHECK(splicing_matrix_resize(samples, noiso, noSamples));
  SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));

  SPLICING_CHECK(splicing_matrix_init(&chainMeans, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &chainMeans);
  SPLICING_CHECK(splicing_matrix_init(&chainVars, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &chainVars);

  /* Initialize Psi(0) randomly */

  SPLICING_CHECK(splicing_drift_proposal_init(noiso, noChains, psi, alpha, &sigma,
					 start, start_psi, gff,
					 gene, readLength, overHang, 
					 position, cigarstr, /*paired=*/ 1, 
					 fragmentProb, fragmentStart,
					 normalMean, normalVar, numDevs));

  SPLICING_CHECK(splicing_drift_proposal_propose(noiso, noChains, 
						 alpha, sigma, psi, alpha));
  
  /* Initialize assignments of reads */  
  
  SPLICING_CHECK(splicing_reassign_samples_paired(mymatch_matrix,
						  &match_order, 
						  psi, noiso, noChains, 
						  fragmentStart,
						  &vass));
  
  /* foreach Iteration m=1, ..., M do */

  while (1) { 

    for (m=0; m < noIterations; m++) {
      
      SPLICING_CHECK(splicing_drift_proposal_propose(noiso, noChains,
						     alpha, sigma,
						     psiNew, alphaNew));

      SPLICING_CHECK(splicing_metropolis_hastings_ratio_paired(&vass,
					     noReads, noChains,
					     psiNew, alphaNew, psi, alpha,
					     sigma, noiso, &isolen, hyperp, 
					     &isoscores, &assscores,
					     &fragmentLength, fragmentStart,
					     m > 0 ? 1 : 0, 
					     &acceptP, &cJS, &pJS));

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
	  memcpy(VECTOR(*logLik)+noS, VECTOR(cJS), noChains * sizeof(double));
	  noS += noChains;
	  lagCounter = 0;
	} else {
	  lagCounter ++;
	}
      }
      
      SPLICING_CHECK(splicing_reassign_samples_paired(mymatch_matrix,
						      &match_order,
						      psi, noiso, noChains, 
						      fragmentStart, &vass));

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
  splicing_matrix_destroy(&valphaNew);
  splicing_matrix_destroy(&valpha);
  splicing_matrix_destroy(&vpsiNew);
  splicing_matrix_destroy(&vpsi);
  splicing_matrix_int_destroy(&vass);
  splicing_vector_destroy(&cJS);
  splicing_vector_destroy(&pJS);
  splicing_vector_destroy(&acceptP);
  SPLICING_FINALLY_CLEAN(8);

  if (!fragmentProb) { 
    splicing_vector_destroy(&vfragmentProb);
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}

int splicing_i_miso_classes1(const splicing_matrix_t *match_matrix,
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

    SPLICING_CHECK(splicing_matrix_resize(class_templates, noiso, 1));
    SPLICING_CHECK(splicing_vector_resize(class_counts, 1));
    
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

#define SETCOL(to, from) do {						\
  int j;								\
  for (j=0; j<noiso; j++) {						\
    MATRIX(*class_templates, j, to) =					\
      MATRIX(*match_matrix, j, from) != 0;				\
  } } while (0)    

int splicing_i_miso_classes2(const splicing_matrix_t *match_matrix,
			     splicing_matrix_t *class_templates,
			     splicing_vector_t *class_counts) {

  int noiso=splicing_matrix_nrow(match_matrix);
  
  if (splicing_matrix_size(match_matrix) == 0) { 

    /* Special case: no reads */
    SPLICING_CHECK(splicing_matrix_resize(class_templates, noiso, 0));
    SPLICING_CHECK(splicing_vector_resize(class_counts, 0));

  } else { 

    int i, j, noreads=splicing_matrix_ncol(match_matrix);
    int lastclass=0;
    double *prev, *curr;
    splicing_vector_int_t match_order;
    int *order;

    SPLICING_CHECK(splicing_vector_int_init(&match_order, noreads));
    SPLICING_FINALLY(splicing_vector_int_destroy, &match_order);

    SPLICING_CHECK(splicing_matrix_binorder_cols(match_matrix, &match_order));
    order=VECTOR(match_order);
 
    SPLICING_CHECK(splicing_matrix_resize(class_templates, noiso, 1));
    SPLICING_CHECK(splicing_vector_resize(class_counts, 1));
   
    SETCOL(lastclass, order[0]);
    prev=&MATRIX(*class_templates, 0, lastclass);
    
    for (i=0; i<noreads; i++) {
      int same=1;
      curr=&MATRIX(*match_matrix, 0, order[i]);
      for (j=0; same && j<noiso; j++) { 
	same=(prev[j] && curr[j]) || (!prev[j] && !curr[j]);
      }
      if (!same) {
	SPLICING_CHECK(splicing_matrix_add_cols(class_templates, 1));
	SPLICING_CHECK(splicing_vector_push_back(class_counts, 0));
	lastclass++;
	prev = &MATRIX(*class_templates, 0, lastclass);
	SETCOL(lastclass, order[i]);
      }
      VECTOR(*class_counts)[lastclass] += 1;
    }
    
    splicing_vector_int_destroy(&match_order);
    SPLICING_FINALLY_CLEAN(1);
  }  
  
  return 0;
}

#undef SETCOL

int splicing_i_miso_classes(const splicing_matrix_t *match_matrix,
			    const splicing_vector_int_t *match_order,
			    splicing_matrix_t *class_templates,
			    splicing_vector_t *class_counts, 
			    splicing_matrix_t *bin_class_templates,
			    splicing_vector_t *bin_class_counts) {

  if (class_templates) { 
    SPLICING_CHECK(splicing_i_miso_classes1(match_matrix, match_order,
					    class_templates, class_counts));
  }
  if (bin_class_templates) { 
    SPLICING_CHECK(splicing_i_miso_classes2(match_matrix,
					    bin_class_templates, 
					    bin_class_counts));
  }
  return 0;
}

