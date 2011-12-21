
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
			    splicing_strvector_t *cigar, 
			    splicing_vector_t *sample_prob) {
  
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

  if (sample_prob) {
    SPLICING_CHECK(splicing_vector_update(sample_prob, &sampleprob));
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

  /* TODO: some of this is already calculated */
  SPLICING_CHECK(splicing_iso_to_genomic(gff, gene, isoform, /*converter=*/ 0,
					 position));

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
  int fragmentEnd;
  int i, j;

  *fragmentStart = normalMean - normalSd * numDevs;
  fragmentEnd    = normalMean + normalSd * numDevs;
  if (*fragmentStart < minLength) { *fragmentStart = minLength; }
  if (fragmentEnd    < *fragmentStart) { fragmentEnd = *fragmentStart; }
  
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
				   splicing_strvector_t *cigar, 
				   splicing_vector_t *sampleprob) {
  
  size_t i, j, noiso, il, nogenes;
  splicing_vector_t *mysampleprob=sampleprob, vsampleprob;
  splicing_vector_t px, cpx;
  double sumpx, sumpsi=0.0;
  splicing_vector_int_t isolen;
  int goodiso=0;
  splicing_vector_int_t exstart, exend, exidx;
  splicing_vector_t *myfragmentProb=(splicing_vector_t*) fragmentProb,
    vfragmentProb;
  int fs, fl;

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
  fs=fragmentStart;
  fl=fragmentStart+il-1;

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));
    
  if ( fabs(splicing_vector_sum(myfragmentProb) - 1.0) > 1e-10 ) {
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

  if (!sampleprob) {
    mysampleprob=&vsampleprob;
    SPLICING_CHECK(splicing_vector_init(mysampleprob, noiso));
    SPLICING_FINALLY(splicing_vector_destroy, mysampleprob);
  } else {
    SPLICING_CHECK(splicing_vector_resize(mysampleprob, noiso));
  }

  for (sumpx=VECTOR(px)[0], i=1; i<il; i++) {
    VECTOR(px)[i] += VECTOR(px)[i-1];
    sumpx += VECTOR(px)[i];
  }
  VECTOR(cpx)[0] = VECTOR(px)[0];
  for (i=1; i<il; i++) {
    VECTOR(cpx)[i] = VECTOR(cpx)[i-1] + VECTOR(px)[i];
  }

  for (i=0; i<noiso; i++) {
    int ilen=VECTOR(isolen)[i];
    int r1= ilen >= fl ? ilen - fl + 1 : 0;
    int r2= ilen >= fs ? (ilen >= fl ? fl - fs : ilen - fs + 1) : 0;
    /* int r3= fs - 1; */
    double sp=0.0;
    if (r1 > 0) { sp += r1; } 
    if (r2 > 0) { sp += VECTOR(cpx)[r2-1]; }
    VECTOR(*mysampleprob)[i] = sp * VECTOR(*expression)[i];
    if (VECTOR(*mysampleprob)[i] != 0) { goodiso += 1; }
    sumpsi += VECTOR(*mysampleprob)[i];
  }

  if (goodiso == 0) {
    SPLICING_ERROR("No isoform is possible", SPLICING_FAILURE);
  }

  for (i=1; i<noiso; i++) {
    VECTOR(*mysampleprob)[i] += VECTOR(*mysampleprob)[i-1];
  }

  SPLICING_CHECK(splicing_vector_int_resize(isoform, noreads*2));

  for (i=0; i<2*noreads; i+=2) {
    int w;
    double rand;
    if (noiso==1) {
      w=0;
    } else if (noiso==2) {
      rand = RNG_UNIF01() * sumpsi;
      w = (rand < VECTOR(*mysampleprob)[0]) ? 0 : 1;
    } else {
      rand = RNG_UNIF01() * sumpsi;
      for (w=0; rand > VECTOR(*mysampleprob)[w]; w++) ;
    }
    VECTOR(*isoform)[i]=VECTOR(*isoform)[i+1]=w;
  }

  if (!sampleprob) { 
    splicing_vector_destroy(mysampleprob);
    SPLICING_FINALLY_CLEAN(1);
  } else {
    for (i=noiso-1; i>0; i--) {
      VECTOR(*mysampleprob)[i] -= VECTOR(*mysampleprob)[i-1];
    }
  }

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
    int ilen=VECTOR(isolen)[iso];
    int r1= ilen >= fl ? ilen - fl + 1 : 0;
    int r2= ilen >= fs ? (ilen >= fl ? fl - fs : ilen - fs + 1) : 0;
    /* int r3= fs - 1; */
    int pos, fragment;
    double sp=0.0;
    if (r1 > 0) { sp += r1; } 
    if (r2 > 0) { sp += VECTOR(cpx)[r2-1]; }
    double rand=RNG_UNIF(0, sp);
    if (rand < r1) { 
      pos = ceil(rand);
    } else {
      int w;
      rand -= r1;
      for (w=0; VECTOR(cpx)[w] < rand; w++) ;
      pos = r1 + r2 - w;
    }

    if (pos <= r1) {
      rand=RNG_UNIF(0, 1.0);
    } else {
      rand=RNG_UNIF(0, VECTOR(px)[r1+r2-pos]);
    }
    for (fragment=0; VECTOR(px)[fragment] < rand; fragment++) ;
    fragment += fragmentStart;

    VECTOR(*position)[j++] = pos;
    VECTOR(*position)[j++] = pos+fragment-readLength;
    
  }

  /* Translate positions to genomic coordinates */

  /* TODO: some of this is already calculated */
  SPLICING_CHECK(splicing_iso_to_genomic(gff, gene, isoform, /*converter=*/ 0,
					 position));

  /* CIGAR strings */

  splicing_strvector_clear(cigar);
  SPLICING_CHECK(splicing_strvector_reserve(cigar, 2*noreads));
  for (j=0; j<2*noreads; j++) {
    char tmp[1000], *tmp2=tmp;
    int iso=VECTOR(*isoform)[j];
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
