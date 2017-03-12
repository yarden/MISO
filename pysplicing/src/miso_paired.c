
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
	sumpsi += MATRIX(*psi, j, k) * (*cptr);				\
	VECTOR(validIso)[noValid] = j;					\
	VECTOR(cumsum)[noValid] = sumpsi;				\
	noValid++;							\
      }									\
    }									\
  } while (0)

int splicing_reassign_samples_paired1(
			     const splicing_matrix_t *matches, 
			     const splicing_vector_int_t *match_order,
			     const splicing_matrix_t *psi, 
			     int noiso, int noChains, int fragmentStart, 
			     splicing_matrix_int_t *result);

int splicing_reassign_samples_paired(
			     const splicing_vector_ptr_t *all_matches, 
			     const splicing_vector_ptr_t *all_match_order,
			     const splicing_vector_ptr_t *all_psi, 
			     int noiso, int noChains, int fragmentStart, 
			     splicing_vector_ptr_t *all_result) {

  int rep, norep = (int) splicing_vector_ptr_size(all_matches);
  for (rep = 0; rep < norep; rep++) {
    const splicing_matrix_t *matches = VECTOR(*all_matches)[rep];
    const splicing_vector_int_t *match_order = VECTOR(*all_match_order)[rep];
    const splicing_matrix_t *psi = VECTOR(*all_psi)[rep];
    splicing_matrix_int_t *result = VECTOR(*all_result)[rep];
    splicing_reassign_samples_paired1(matches, match_order, psi, noiso, noChains,
				      fragmentStart, result);
  }

  return 0;
}  

int splicing_reassign_samples_paired1(
			     const splicing_matrix_t *matches, 
			     const splicing_vector_int_t *match_order,
			     const splicing_matrix_t *psi, 
			     int noiso, int noChains, int fragmentStart, 
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

int splicing_score_iso_paired(const splicing_vector_t *psi, int noiso, 
			      const splicing_vector_int_t *assignment, 
			      const splicing_vector_int_t *pisolen, 
			      const splicing_vector_t *assscores,
			      double *res) {

  int noreads = (int) splicing_vector_int_size(assignment);
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
				splicing_miso_hyperprior_t *hyperprior,
				const splicing_vector_int_t *isolen,
				const splicing_matrix_t *isoscores, 
				const splicing_vector_t *assscores,
				const splicing_matrix_int_t *fragmentLength,
				int fragmentStart, 
				splicing_vector_t *score) {

  int i, j, noiso = (int) splicing_vector_int_size(isolen);

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
    if (hyperprior->prior == SPLICING_MISO_PRIOR_DIRICHLET) {
      SPLICING_CHECK(splicing_ldirichlet(&tmp, &(hyperprior->dirichlet_hyperp),
					 noiso, &psiProb));
    } else {
      SPLICING_CHECK(splicing_llogistic(&tmp,
					VECTOR(hyperprior->logistic_mean)[j],
					VECTOR(hyperprior->logistic_var)[j],
					&psiProb));
    }
    
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
			      splicing_miso_hyperprior_t *hyperprior,
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
					     psiNew, hyperprior,
					     isolen, isoscores, assscores, 
					     fragmentLength, fragmentStart,
					     ppJS));
  SPLICING_CHECK(splicing_score_joint_paired(ass, no_reads, noChains, 
					     psi, hyperprior,
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
			 const splicing_replicate_reads_t *reads,
			 int readLength, int overHang,
			 int noChains, int noIterations, 
			 int maxIterations, int noBurnIn, int noLag,
			 splicing_miso_hyperprior_t *hyperprior,
			 splicing_miso_start_t start, 
			 splicing_miso_stop_t stop,
			 const splicing_matrix_t *start_psi,
			 const splicing_vector_t *fragmentProb,
			 int fragmentStart, double normalMean, 
			 double normalVar, double numDevs,

			 splicing_matrix_t *pop_samples,
			 splicing_vector_ptr_t *samples, /* matrix_t */
			 
			 splicing_vector_t *logLik,
			 splicing_vector_ptr_t *match_matrix, /* matrix_t */
			 splicing_matrix_t *class_templates,
			 splicing_vector_ptr_t *class_counts, /* vector_t */
			 splicing_matrix_t *bin_class_templates,
			 splicing_vector_ptr_t *bin_class_counts, /* vector_t */
			 splicing_vector_ptr_t *assignment,	  /* vector_int_t */
			 splicing_miso_rundata_t *rundata) {

  size_t noiso;
  splicing_vector_t acceptP, cJS, pJS;
  double sigma;
  splicing_vector_int_t noReads;
  int noReplicates = splicing_replicate_reads_noreps(reads);
  splicing_vector_ptr_t vass;	/* matrix_int_t */

  /* These all contain matrix_t's */
  splicing_vector_ptr_t vpsi, vpsiNew, valpha, valphaNew, 
    *psi=&vpsi, *psiNew=&vpsiNew, *alpha=&valpha, *alphaNew=&valphaNew;

  int noSamples = noChains * (noIterations - noBurnIn) / noLag;
  int i, j, m, lagCounter=0, noS=0;
  splicing_vector_ptr_t *mymatch_matrix=match_matrix, vmatch_matrix; /* matrix_t */
  splicing_vector_ptr_t match_order;				     /* vector_int_t */
  splicing_vector_int_t isolen;
  splicing_matrix_t isoscores;
  splicing_vector_ptr_t fragmentLength; /* matrix_int_t */
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
					    readLength, myfragmentProb,
					    &fragmentStart));
    splicing_vector_scale(myfragmentProb, 
			  1.0/splicing_vector_sum(myfragmentProb));
  }

  il = (int) splicing_vector_size(myfragmentProb);

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

  if (splicing_vector_min(&hyperprior->logistic_var) < 0) {
    SPLICING_ERROR("Variance of the logistic prior must be non-negative",
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

  if (!pop_samples && !samples) {
    SPLICING_ERROR("At least one of 'samples' and 'pop_samples' must"
		   "be non-NULL", SPLICING_EINVAL);
  }

  if (noReplicates > 1 && ! pop_samples) {
    SPLICING_ERROR("Replicates, but pop_samples is NULL", SPLICING_EINVAL);
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
    VECTOR(noReads)[i] = splicing_replicate_reads_noreads(reads, i) / 2;
  }

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

  SPLICING_CHECK(splicing_vector_ptr_init(&fragmentLength, noReplicates));
  splicing_vector_ptr_set_item_destructor(&fragmentLength,
    (splicing_finally_func_t *) splicing_matrix_int_destroy_free);
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &fragmentLength);
  for (i = 0; i < noReplicates; i++) {
    splicing_matrix_int_t *fl = splicing_Calloc(1, splicing_matrix_int_t);
    if (!fl) { SPLICING_ERROR("No memory for fragment length", SPLICING_ENOMEM); }
    VECTOR(fragmentLength)[i] = fl;
    SPLICING_CHECK(splicing_matrix_int_init(fl, noiso, VECTOR(noReads)[i]));
  }

  for (i = 0; i < noReplicates; i++) {
    const splicing_vector_int_t *position =
      splicing_replicate_reads_pos(reads, i);
    const char **cigarstr = splicing_replicate_reads_cigar(reads, i);
    splicing_matrix_t *mm = VECTOR(*mymatch_matrix)[i];
    splicing_vector_int_t *mo = VECTOR(match_order)[i];
    splicing_matrix_int_t *fl = VECTOR(fragmentLength)[i];
    SPLICING_CHECK(splicing_matchIso_paired(gff, (int) gene, position, cigarstr,
					    readLength, overHang, 
					    myfragmentProb, 
					    fragmentStart, normalMean, 
					    normalVar, numDevs, mm, fl));
    SPLICING_CHECK(splicing_order_matches(mm, mo));
  }
  
  if (class_templates || bin_class_templates) {
    if (class_counts) {
      SPLICING_CHECK(splicing_vector_ptr_resize(class_counts, noReplicates));
      splicing_vector_ptr_set_item_destructor(class_counts,
        (splicing_finally_func_t *) splicing_vector_destroy_free);
      SPLICING_FINALLY(splicing_vector_ptr_destroy, class_counts);
      for (i = 0; i < noReplicates; i++) {
	splicing_vector_t *cc = splicing_Calloc(1, splicing_vector_t);
	if (!cc) { SPLICING_ERROR("No memory for class counts", SPLICING_ENOMEM); }
	VECTOR(*class_counts)[i] = cc;
	SPLICING_CHECK(splicing_vector_init(cc, 0));
      }
    }
    if (bin_class_counts) {
      SPLICING_CHECK(splicing_vector_ptr_resize(bin_class_counts, noReplicates));
      splicing_vector_ptr_set_item_destructor(bin_class_counts,
        (splicing_finally_func_t *) splicing_vector_destroy_free);
      SPLICING_FINALLY(splicing_vector_ptr_destroy, bin_class_counts);
      for (i = 0; i < noReplicates; i++) {
	splicing_vector_t *bcc = splicing_Calloc(1, splicing_vector_t);
	if (!bcc) { SPLICING_ERROR("No memory for class counts", SPLICING_ENOMEM); }
	VECTOR(*bin_class_counts)[i] = bcc;
	SPLICING_CHECK(splicing_vector_init(bcc, 0));
      }
    }

    for (i = 0; i < noReplicates; i++) {
      splicing_matrix_t *mm = VECTOR(*mymatch_matrix)[i];
      splicing_vector_int_t *mo = VECTOR(match_order)[i];
      splicing_vector_t *cc = class_counts ? VECTOR(*class_counts)[i] : 0;
      splicing_vector_t *bcc = bin_class_counts ? VECTOR(*bin_class_counts)[i] : 0;
      SPLICING_CHECK(splicing_i_miso_classes(mm, mo, class_templates, cc,
					     bin_class_templates,  bcc));
    }
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

  SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));

  SPLICING_CHECK(splicing_matrix_init(&chainMeans, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &chainMeans);
  SPLICING_CHECK(splicing_matrix_init(&chainVars, noiso, noChains));
  SPLICING_FINALLY(splicing_matrix_destroy, &chainVars);

  if (pop_samples) {
    SPLICING_CHECK(splicing_matrix_resize(pop_samples, noiso, noSamples));
  }

  if (samples) {
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
  }  
  
  /* Initialize Psi(0) randomly */

  SPLICING_CHECK(splicing_drift_proposal_init_paired((int)noiso, noChains, psi, alpha,
					 &sigma, start, start_psi, gff,
					 (int) gene, readLength, overHang,
					 reads, /*paired=*/ 1, 
					 fragmentProb, fragmentStart,
					 normalMean, normalVar, numDevs));

  SPLICING_CHECK(splicing_drift_proposal_propose_paired((int) noiso, noChains,
						 alpha, sigma, psi, alpha));
  
  /* Initialize assignments of reads */  
  
  SPLICING_CHECK(splicing_reassign_samples_paired(mymatch_matrix,
						  &match_order, 
						  psi, (int) noiso, noChains,
						  fragmentStart,
						  &vass));
  
  /* foreach Iteration m=1, ..., M do */

  while (1) { 

    for (m=0; m < noIterations; m++) {
      
      SPLICING_CHECK(splicing_drift_proposal_propose_paired((int) noiso, noChains,
						     alpha, sigma,
						     psiNew, alphaNew));

      for (i = 0; i < noReplicates; i++) {

	splicing_matrix_int_t *vass1 = VECTOR(vass)[i];
	splicing_matrix_t *psiNew1 = VECTOR(*psiNew)[i];
	splicing_matrix_t *alphaNew1 = VECTOR(*alphaNew)[i];
	splicing_matrix_t *psi1 = VECTOR(*psi)[i];
	splicing_matrix_t *alpha1 = VECTOR(*alpha)[i];
	splicing_matrix_int_t *fragmentLength1 = VECTOR(fragmentLength)[i];
	

	SPLICING_CHECK(splicing_metropolis_hastings_ratio_paired(vass1,
 	  VECTOR(noReads)[i], noChains, psiNew1, alphaNew1, psi1, alpha1,
	  sigma, (int) noiso, &isolen, hyperprior,
	  &isoscores, &assscores, fragmentLength1, fragmentStart,
	  m > 0 ? 1 : 0, &acceptP, &cJS, &pJS));

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

      }	/*  i < noReplicates */

      /* Update population mean, var */
      if (noReplicates > 1) {
	splicing_update_replicate_mean(hyperprior, psi);
	splicing_update_replicate_var(hyperprior, psi);
      }
      
      if (m >= noBurnIn) {
	if (lagCounter == noLag - 1) {

	  memcpy(VECTOR(*logLik)+noS, VECTOR(cJS), noChains * sizeof(double));
	  
	  if (samples) {
	    for (i = 0; i < noReplicates; i++) {
	      splicing_matrix_t *samples1 = VECTOR(*samples)[i];
	      splicing_matrix_t *psi1 = VECTOR(*psi)[i];
	      memcpy(&MATRIX(*samples1, 0, noS), &MATRIX(*psi1, 0, 0), 
		     noChains * noiso * sizeof(double));
	    }
	  }

	  if (pop_samples) {
	    if (noReplicates == 1 && samples) {
	      splicing_matrix_t *psi1 = VECTOR(*psi)[0];
	      memcpy(&MATRIX(*pop_samples, 0, noS), &MATRIX(*psi1, 0, 0),
		     noChains * noiso * sizeof(double));
	    } else {
	      splicing_logit_inv_raw(VECTOR(hyperprior->logistic_mean),
				     &MATRIX(*pop_samples, 0, noS),
				     (int) noiso - 1, noChains);
	    }
	  }

	  noS += noChains;
	  lagCounter = 0;
	  
	} else {
	  lagCounter ++;
	}
      }
      
      SPLICING_CHECK(splicing_reassign_samples_paired(mymatch_matrix,
						      &match_order,
						      psi, (int) noiso, noChains,
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
	splicing_matrix_t *sam = pop_samples ? pop_samples : VECTOR(*samples)[0];
	SPLICING_CHECK(splicing_i_check_convergent_mean(&chainMeans,
							&chainVars, sam,
							&shouldstop));
      }
      break;
    }

    if (shouldstop) { break; }
    
    noS=0;
    noIterations = 3*noIterations - 2*noBurnIn;
    noBurnIn = m;
    noSamples = noChains * (noIterations - noBurnIn) / noLag;
    lagCounter = 0;

    if (samples) {
      for (i = 0; i < noReplicates; i++) {
	splicing_matrix_t *sam = VECTOR(*samples)[i];
	SPLICING_CHECK(splicing_matrix_resize(sam, noiso, noSamples));
      }
    }
    if (pop_samples) {
      SPLICING_CHECK(splicing_matrix_resize(pop_samples, noiso, noSamples));
    }
    SPLICING_CHECK(splicing_vector_resize(logLik, noSamples));
  }

  splicing_matrix_destroy(&chainVars);
  splicing_matrix_destroy(&chainMeans);
  SPLICING_FINALLY_CLEAN(2);

  if (assignment) {
    splicing_vector_ptr_resize(assignment, noReplicates);
    splicing_vector_ptr_set_item_destructor(assignment,
      (splicing_finally_func_t *) splicing_vector_int_destroy_free);
    SPLICING_FINALLY(splicing_vector_ptr_destroy, assignment);

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

  splicing_vector_destroy(&assscores);
  splicing_matrix_destroy(&isoscores);
  splicing_vector_int_destroy(&isolen);
  splicing_vector_ptr_destroy(&fragmentLength);
  splicing_vector_ptr_destroy(&match_order);
  SPLICING_FINALLY_CLEAN(5);
  if (!match_matrix) {
    splicing_vector_ptr_destroy(mymatch_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  splicing_vector_ptr_destroy(&valphaNew);
  splicing_vector_ptr_destroy(&valpha);
  splicing_vector_ptr_destroy(&vpsiNew);
  splicing_vector_ptr_destroy(&vpsi);
  splicing_vector_ptr_destroy(&vass);
  splicing_vector_destroy(&cJS);
  splicing_vector_destroy(&pJS);
  splicing_vector_destroy(&acceptP);
  SPLICING_FINALLY_CLEAN(8);

  if (!fragmentProb) { 
    splicing_vector_destroy(&vfragmentProb);
    SPLICING_FINALLY_CLEAN(1);
  }

  /* We always return the same number of samples. */

  if (rundata->noSamples != noSamples) {
    splicing_vector_remove_section(logLik, 0, noSamples-rundata->noSamples);
    if (samples) {
      for (i = 0; i < noReplicates; i++) {
	splicing_matrix_t *sam = VECTOR(*samples)[i];
	splicing_matrix_remove_cols_section(sam, 0,
					  noSamples-rundata->noSamples);
      }
    }
    if (pop_samples) {
      splicing_matrix_remove_cols_section(pop_samples, 0,
					  noSamples-rundata->noSamples);
    }
  }
  
  return 0;
}

int splicing_i_miso_classes1(const splicing_matrix_t *match_matrix,
			     const splicing_vector_int_t *match_order,
			     splicing_matrix_t *class_templates,
			     splicing_vector_t *class_counts) {
  
  int noiso = (int) splicing_matrix_nrow(match_matrix);
  
  if (splicing_matrix_size(match_matrix) == 0) { 

    /* Special case: no reads */
    splicing_matrix_resize(class_templates, noiso, 0);
    splicing_vector_resize(class_counts, 0);

  } else { 

    int i, noreads = (int) splicing_vector_int_size(match_order);
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

  int noiso = (int) splicing_matrix_nrow(match_matrix);
  
  if (splicing_matrix_size(match_matrix) == 0) { 

    /* Special case: no reads */
    SPLICING_CHECK(splicing_matrix_resize(class_templates, noiso, 0));
    SPLICING_CHECK(splicing_vector_resize(class_counts, 0));

  } else { 

    int i, j, noreads = (int) splicing_matrix_ncol(match_matrix);
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

