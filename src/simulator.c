
#include "splicing.h"
#include "splicing_random.h"
#include "splicing_error.h"

#include <stdio.h>
#include <math.h>

int splicing_create_gene(const splicing_vector_int_t *exons,
			 const splicing_vector_int_t *isoforms,
			 const char *id, const char *seqid, 
			 const char *source, splicing_strand_t strand,
			 splicing_gff_t *extend) {

  size_t i=0;
  size_t exlen=splicing_vector_int_size(exons);
  size_t isolen=splicing_vector_int_size(isoforms);
  size_t genestart=splicing_vector_int_min(exons);
  size_t geneend=splicing_vector_int_max(exons);
  char buffer[5000], buffer2[5000];
  int noiso=0;
  
  /* TODO: error checks */
  
  /* Gene */
  SPLICING_CHECK(splicing_gff_append(extend, seqid, source, 
				     SPLICING_TYPE_GENE, 
				     genestart, geneend, 
				     /*score=*/ SPLICING_NA_REAL, 
				     strand, /*phase=*/ SPLICING_NA_INTEGER,
				     id, /*parent=*/ 0));

  while (i<isolen) {
    size_t mmin=VECTOR(*exons)[ 2*VECTOR(*isoforms)[i] ];
    size_t mmax=VECTOR(*exons)[ 2*VECTOR(*isoforms)[i] + 1 ];
    size_t j, exon=0;    
    for (j=i+1; VECTOR(*isoforms)[j] >= 0; j++) {
      size_t m1=VECTOR(*exons)[ 2*VECTOR(*isoforms)[j] ];
      size_t m2=VECTOR(*exons)[ 2*VECTOR(*isoforms)[j] + 1 ];      
      if (m1 < mmin) { mmin = m1; }
      if (m2 > mmax) { mmax = m2; }
    }
    snprintf(buffer, sizeof(buffer)/sizeof(char)-sizeof(char), 
	     "%s-isoform-%i", id, noiso);    
    SPLICING_CHECK(splicing_gff_append(extend, seqid, source, 
				       SPLICING_TYPE_MRNA, mmin, mmax, 
				       /*score=*/ SPLICING_NA_REAL, strand,
				       /*phase=*/ SPLICING_NA_INTEGER,
				       buffer, /*parent=*/ id));
    for (; VECTOR(*isoforms)[i] >= 0; i++) {
      snprintf(buffer2, sizeof(buffer2)/sizeof(char)-sizeof(char),
	       "%s-isoform-%i-exon-%i", id, (int) noiso, (int) exon++);
      SPLICING_CHECK(splicing_gff_append(extend, seqid, source, 
			 SPLICING_TYPE_EXON,
			 VECTOR(*exons)[ 2*VECTOR(*isoforms)[i] ],
			 VECTOR(*exons)[ 2*VECTOR(*isoforms)[i] + 1 ],
			 /*score=*/ SPLICING_NA_REAL, strand,
			 /*phase=*/ SPLICING_NA_INTEGER, buffer2, 
			 /*parent=*/ buffer));
    }
    noiso++;
    i++;
  }
  
  return 0;
}

int splicing_simulate_reads(const splicing_gff_t *gff, int gene,
			    const splicing_vector_t *expression,
			    int noreads, int readLength,
			    splicing_vector_int_t *isoform, 
			    splicing_vector_int_t *position, 
			    splicing_strvector_t *cigar) {
  
  size_t i, p, noiso, goodiso=0, nogenes;
  splicing_vector_int_t effisolen;
  splicing_vector_t sampleprob;
  double rand, sumpsi=0.0;
  splicing_vector_int_t exstart, exend, exidx;

  SPLICING_CHECK(splicing_gff_nogenes(gff, &nogenes));
  if (gene < 0 || gene >= nogenes) {
    SPLICING_ERROR("Invalid gene id", SPLICING_EINVAL);
  }

  /* TODO: more error checks */

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));
    
  SPLICING_CHECK(splicing_vector_int_init(&effisolen, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &effisolen);
  SPLICING_CHECK(splicing_vector_init(&sampleprob, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &sampleprob);
  SPLICING_CHECK(splicing_vector_int_resize(isoform, noreads));
  SPLICING_CHECK(splicing_gff_isolength_one(gff, gene, &effisolen));
  for (i=0; i<noiso; i++) {
    int l=VECTOR(effisolen)[i]-readLength+1;
    VECTOR(effisolen)[i] = l > 0 ? l : 0;
    VECTOR(sampleprob)[i] = VECTOR(*expression)[i] * VECTOR(effisolen)[i];
    if (VECTOR(sampleprob)[i] != 0) { goodiso++; }
    sumpsi += VECTOR(sampleprob)[i];
  }

  if (goodiso==0) {
    SPLICING_ERROR("No isoform is possible", SPLICING_FAILURE);
  }

  for (i=1; i<noiso; i++) {
    VECTOR(sampleprob)[i] += VECTOR(sampleprob)[i-1];
  }

  for (i=0; i<noreads; i++) {
    int w;
    if (noiso==1) {
      w=0;
    } else if (noiso==2) {
      rand = RNG_UNIF01() * sumpsi;
      w = (rand < VECTOR(sampleprob)[0]) ? 0 : 1;
    } else {
      rand = RNG_UNIF01() * sumpsi;
      for (w=0; rand > VECTOR(sampleprob)[w]; w++) ;
    }
    VECTOR(*isoform)[i]=w;
  }
  
  splicing_vector_destroy(&sampleprob);
  SPLICING_FINALLY_CLEAN(1);

  /* OK, we have the isoforms, now we need the read positions, 
     these are uniformly sampled from the individual isoforms. */

  SPLICING_CHECK(splicing_vector_int_resize(position, noreads));
  SPLICING_CHECK(splicing_vector_int_init(&exstart, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exstart);
  SPLICING_CHECK(splicing_vector_int_init(&exend, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exend);
  SPLICING_CHECK(splicing_vector_int_init(&exidx, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exidx);
  SPLICING_CHECK(splicing_gff_exon_start_end(gff, &exstart, &exend, &exidx,
					     gene));

  /* Positions in isoform coordinates first */

  for (i=0; i<noreads; i++) { 
    int iso=VECTOR(*isoform)[i];
    int len=VECTOR(effisolen)[iso];
    VECTOR(*position)[i]=RNG_INTEGER(1, len);
  }

  /* Translate isoform coordinates to genomic coordintes */

  SPLICING_CHECK(splicing_iso_to_genomic(gff, gene, isoform, &exstart, 
					 &exend, &exidx, position));

  /* CIGAR strings */

  splicing_strvector_clear(cigar);
  SPLICING_CHECK(splicing_strvector_reserve(cigar, noreads));
  for (i=0; i<noreads; i++) {
    char tmp[1000], *tmp2=tmp;
    int iso=VECTOR(*isoform)[i];
    size_t rs=VECTOR(*position)[i];
    int ex=0;
    int rl=readLength;
    for (ex=VECTOR(exidx)[iso]; VECTOR(exend)[ex] < rs; ex++) ;
    while (VECTOR(exend)[ex] < rs+rl-1) {
      tmp2 += snprintf(tmp2, sizeof(tmp)/sizeof(char)-(tmp2-tmp)-1, "%iM%iN",
		       (int) (VECTOR(exend)[ex]-rs+1), 
		       (int) (VECTOR(exstart)[ex+1]-VECTOR(exend)[ex]-1));
      if (tmp2 >= tmp + sizeof(tmp)/sizeof(char)) {
	SPLICING_ERROR("CIGAR string too long", SPLICING_EINVAL);
      }
      rl -= (VECTOR(exend)[ex] - rs + 1);
      rs = VECTOR(exstart)[ex+1];
      ex++;
    }
    tmp2 += snprintf(tmp2, sizeof(tmp)/sizeof(char)-(tmp2-tmp)-1, "%iM", rl);
    if (tmp2 >= tmp + sizeof(tmp)/sizeof(char)) {
      SPLICING_ERROR("CIGAR string too long", SPLICING_EINVAL); }
    SPLICING_CHECK(splicing_strvector_append(cigar, tmp));
  }

  splicing_vector_int_destroy(&exidx);
  splicing_vector_int_destroy(&exend);
  splicing_vector_int_destroy(&exstart);
  splicing_vector_int_destroy(&effisolen);
  SPLICING_FINALLY_CLEAN(4);
  
  return 0;
}

int splicing_normal_fragment(double normalMean, double normalVar, 
			     double numDevs, int minLength, 
			     splicing_vector_t *fragmentProb,
			     int *fragmentStart) {

  double normalSd=sqrt(normalVar);
  int fragmentEnd = normalMean + normalSd * numDevs;
  int i, j;

  *fragmentStart = normalMean - normalSd * numDevs;
  if (*fragmentStart < minLength) *fragmentStart = minLength;
  
  SPLICING_CHECK(splicing_vector_resize(fragmentProb,
					fragmentEnd - *fragmentStart + 1));
  for (i=*fragmentStart, j=0; i<=fragmentEnd; i++, j++) {
    VECTOR(*fragmentProb)[j] = splicing_dnorm(i, normalMean, normalSd);
  }
  
  return 0;
}

int splicing_simulate_paired_reads(const splicing_gff_t *gff, int gene,
				   const splicing_vector_t *expression,
				   int noreads, int readLength,
				   const splicing_vector_t *fragmentProb,
				   int fragmentStart, double normalMean,
				   double normalVar, double numDevs,
				   splicing_vector_int_t *isoform,
				   splicing_vector_int_t *position,
				   splicing_strvector_t *cigar) {
  
  size_t i, j, noiso, il, nogenes;
  splicing_vector_t sampleprob;
  splicing_vector_t px, cpx;
  double sumpx, sumpsi=0.0;
  splicing_vector_int_t isolen;
  int goodiso=0;
  splicing_vector_int_t exstart, exend, exidx;
  splicing_vector_t *myfragmentProb=(splicing_vector_t*) fragmentProb,
    vfragmentProb;

  SPLICING_CHECK(splicing_gff_nogenes(gff, &nogenes));
  if (gene < 0 || gene >= nogenes) {
    SPLICING_ERROR("Invalid gene id", SPLICING_EINVAL);
  }

  /* TODO: more error checks */

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
    
  if ( fabs(splicing_vector_sum(myfragmentProb) - 1.0) > 1e-13 ) {
    SPLICING_ERROR("Fragment length distribution does not sum up to 1", 
		   SPLICING_EINVAL);
  }

  SPLICING_CHECK(splicing_vector_int_init(&isolen, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isolen);
  SPLICING_CHECK(splicing_gff_isolength_one(gff, gene, &isolen));
  
  SPLICING_CHECK(splicing_vector_copy(&px, myfragmentProb));
  SPLICING_FINALLY(splicing_vector_destroy, &px);
  SPLICING_CHECK(splicing_vector_init(&cpx, il));
  SPLICING_FINALLY(splicing_vector_destroy, &cpx);
  SPLICING_CHECK(splicing_vector_init(&sampleprob, noiso));
  SPLICING_FINALLY(splicing_vector_destroy, &sampleprob);

  for (sumpx=VECTOR(px)[0], i=1; i<il; i++) {
    VECTOR(px)[i] += VECTOR(px)[i-1];
    sumpx += VECTOR(px)[i];
  }
  VECTOR(cpx)[0] = VECTOR(px)[il-1];
  for (i=1; i<il; i++) {
    VECTOR(cpx)[i] = VECTOR(cpx)[i-1] + VECTOR(px)[il-1-i];
  }

  for (i=0; i<noiso; i++) {
    double sp=0.0;
    int j, efflen=VECTOR(isolen)[i] - 2 * readLength;
    if (il+fragmentStart-2*readLength < efflen) {
      sp=efflen-il-(fragmentStart-2*readLength)+1 + sumpx;
    } else if (efflen != 0) {
      sp=VECTOR(cpx)[efflen-(fragmentStart-2*readLength)];
    }
    VECTOR(sampleprob)[i] = sp * VECTOR(*expression)[i];
    if (VECTOR(sampleprob)[i] != 0) { goodiso += 1; }
    sumpsi += VECTOR(sampleprob)[i];
  }

  if (goodiso == 0) {
    SPLICING_ERROR("No isoform is possible", SPLICING_FAILURE);
  }

  for (i=1; i<noiso; i++) {
    VECTOR(sampleprob)[i] += VECTOR(sampleprob)[i-1];
  }

  SPLICING_CHECK(splicing_vector_int_resize(isoform, noreads*2));

  for (i=0; i<2*noreads; i+=2) {
    int w;
    double rand;
    if (noiso==1) {
      w=0;
    } else if (noiso==2) {
      rand = RNG_UNIF01() * sumpsi;
      w = (rand < VECTOR(sampleprob)[0]) ? 0 : 1;
    } else {
      rand = RNG_UNIF01() * sumpsi;
      for (w=0; rand > VECTOR(sampleprob)[w]; w++) ;
    }
    VECTOR(*isoform)[i]=VECTOR(*isoform)[i+1]=w;
  }

  splicing_vector_destroy(&sampleprob);
  SPLICING_FINALLY_CLEAN(1);

  /* We have the isoforms, now get the read positions. */
  
  SPLICING_CHECK(splicing_vector_int_resize(position, noreads*2));
  SPLICING_CHECK(splicing_vector_int_init(&exstart, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exstart);
  SPLICING_CHECK(splicing_vector_int_init(&exend, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exend);
  SPLICING_CHECK(splicing_vector_int_init(&exidx, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exidx);
  SPLICING_CHECK(splicing_gff_exon_start_end(gff, &exstart, &exend, &exidx,
					     gene));
  
  /* Positions in isoform coordinates first. 
     These are sampled based on the fragment length distribution. */

  for (i=0, j=0; i<noreads; i++) {
    int iso=VECTOR(*isoform)[2*i];
    int efflen=VECTOR(isolen)[iso] - 2 * readLength+1;
    int nounr=efflen-il-(fragmentStart-2*readLength);
    int pos, fragment;
    double rand;
    if (il+(fragmentStart-2*readLength) < efflen) {
      rand=RNG_UNIF(0, sumpx*nounr + sumpx);
      if (rand < sumpx*nounr) {
	pos=RNG_INTEGER(0, nounr-1);
      } else {
	rand -= sumpx*nounr;
	for (pos=0; rand > VECTOR(cpx)[pos]; pos++) ;
	pos += nounr;
      }
    } else {
      rand=RNG_UNIF(0, VECTOR(cpx)[il-1] - 
		    VECTOR(cpx)[il+(fragmentStart-2*readLength)-efflen]);
      rand += VECTOR(cpx)[il+(fragmentStart-2*readLength)-efflen-1];
      for (pos=il+(fragmentStart-2*readLength)-efflen; 
	   rand > VECTOR(cpx)[pos]; pos++) ;
    }
    
    if (pos+il < efflen) { 
      rand=RNG_UNIF(0, 1.0);
    } else {
      rand=RNG_UNIF(0, VECTOR(px)[efflen-pos]);
    }
    for (fragment=0; rand > VECTOR(px)[fragment]; fragment++) ;
    fragment += (fragmentStart-2*readLength);

    VECTOR(*position)[j++] = pos+1;
    VECTOR(*position)[j++] = pos+1+fragment+readLength;
  }
  
  /* Translate positions to genomic coordinates */

  SPLICING_CHECK(splicing_iso_to_genomic(gff, gene, isoform, &exstart,
					 &exend, &exidx, position));

  /* CIGAR strings */

  splicing_strvector_clear(cigar);
  SPLICING_CHECK(splicing_strvector_reserve(cigar, 2*noreads));
  for (j=0; j<2*noreads; j++) {
    char tmp[1000], *tmp2=tmp;
    int iso=VECTOR(*isoform)[j/2];
    size_t rs=VECTOR(*position)[j];
    int ex=0;
    int rl=readLength;
    for (ex=VECTOR(exidx)[iso]; VECTOR(exend)[ex] < rs; ex++) ;
    while (rs + rl - 1 > VECTOR(exend)[ex]) {
      tmp2 += snprintf(tmp2, sizeof(tmp)/sizeof(char)-(tmp2-tmp)-1, "%iM%iN",
		       (int) (VECTOR(exend)[ex]-rs+1), 
		       (int) (VECTOR(exstart)[ex+1]-VECTOR(exend)[ex]-1));
      if (tmp2 >= tmp + sizeof(tmp)/sizeof(char)) {
	SPLICING_ERROR("CIGAR string too long", SPLICING_EINVAL);
      }
      rl -= (VECTOR(exend)[ex] - rs + 1);
      rs = VECTOR(exstart)[ex+1];
      ex++;
    }
    tmp2 += snprintf(tmp2, sizeof(tmp)/sizeof(char)-(tmp2-tmp)-1, "%iM", rl);
    if (tmp2 >= tmp + sizeof(tmp)/sizeof(char)) {
      SPLICING_ERROR("CIGAR string too long", SPLICING_EINVAL);
    }
    SPLICING_CHECK(splicing_strvector_append(cigar, tmp));
  }

  splicing_vector_int_destroy(&exidx);
  splicing_vector_int_destroy(&exend);
  splicing_vector_int_destroy(&exstart);
  splicing_vector_destroy(&cpx);
  splicing_vector_destroy(&px);
  splicing_vector_int_destroy(&isolen);
  SPLICING_FINALLY_CLEAN(6);

  if (!fragmentProb) { 
    splicing_vector_destroy(myfragmentProb); 
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}
