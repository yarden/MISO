
#include "splicing.h"
#include "splicing_vector.h"
#include "splicing_error.h"

#include <ctype.h>

int splicing_matchIso(const splicing_gff_t *gff, int gene, 
		      const splicing_vector_int_t *position, 
		      const char **cigarstr, int overHang, int readLength,
		      splicing_matrix_t *result) {

  int noreads=splicing_vector_int_size(position);
  int r, i;
  splicing_vector_int_t exstart, exend, exidx, cigar, cigaridx, cigarlength;
  size_t noiso;

  if (overHang==0) { overHang=1; }
  if (overHang < 1) {
    SPLICING_ERROR("Overhang length invalid. Must be positive", 
		   SPLICING_EINVAL);
  }
  if (readLength < 0) { 
    SPLICING_ERROR("Read length cannot be negative", SPLICING_EINVAL);
  }

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));
  SPLICING_CHECK(splicing_vector_int_init(&exstart, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exstart);
  SPLICING_CHECK(splicing_vector_int_init(&exend, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exend);
  SPLICING_CHECK(splicing_vector_int_init(&exidx, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exidx);
  SPLICING_CHECK(splicing_gff_exon_start_end(gff, &exstart, &exend, 
					     &exidx, gene));

  SPLICING_CHECK(splicing_vector_int_init(&cigar, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &cigar);
  SPLICING_CHECK(splicing_vector_int_init(&cigaridx, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &cigaridx);
  SPLICING_CHECK(splicing_vector_int_init(&cigarlength, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &cigarlength);
  SPLICING_CHECK(splicing_parse_cigar(cigarstr, noreads, &cigar, &cigaridx, 
				      &cigarlength, readLength));

  SPLICING_CHECK(splicing_matrix_resize(result, noiso, noreads));

  for (r=0; r<noreads; r++) {
    int *mycig=VECTOR(cigar) + VECTOR(cigaridx)[r];
    int nocig=VECTOR(cigaridx)[r+1] - VECTOR(cigaridx)[r];
    int len =VECTOR(cigarlength)[r];

    /* If the read length is specified, then filter out the reads
       that are shorter */
    if (len < readLength) { 
      for (i=0; i<noiso; i++) { MATRIX(*result, i, r)=0; }

    /* We can also filter out the reads that do not satisfy the overhang 
       constraint here. We assume that the CIGAR string starts and ends 
       with a match (of non-zero length). */
    } else if (mycig[0] < overHang || mycig[nocig-1] < overHang) {
      for (i=0; i<noiso; i++) { MATRIX(*result, i, r)=0; }
    } else {

      for (i=0; i<noiso; i++) {
	int c, pos=VECTOR(*position)[r];
	int ex=VECTOR(exidx)[i];
	
	/* Look for the exon where the read starts */
	while (ex < VECTOR(exidx)[i+1] &&
	       (pos < VECTOR(exstart)[ex] || VECTOR(exend)[ex] < pos)) {
	  ex++;
	}
	if (ex >= VECTOR(exidx)[i+1]) { MATRIX(*result, i, r)=0; continue; }
	
	/* Got it, match cigar string to exons */
	MATRIX(*result, i, r)=1;
	for (c=0; c<nocig; c++) {
	  if (mycig[c] > 0) { /* exon */
	    if (pos + mycig[c] - 1 > VECTOR(exend)[ex]) { 
	      MATRIX(*result, i, r)=0; break;
	    }
	    pos += mycig[c];
	  } else {	  	   /* intron */
	    if (pos != VECTOR(exend)[ex]+1) {
	      MATRIX(*result, i, r)=0; break; 
	    }
	    pos -= mycig[c];
	    ex += 1;
	    if (ex >= VECTOR(exidx)[i+1] || pos != VECTOR(exstart)[ex]) {
	      MATRIX(*result, i, r)=0; break;
	    }
	  }
	}
      } /* i < noiso */
    }   /* r < noreads */
  }	/* if overhang os OK */

  splicing_vector_int_destroy(&cigarlength);
  splicing_vector_int_destroy(&cigaridx);
  splicing_vector_int_destroy(&cigar);
  splicing_vector_int_destroy(&exidx);
  splicing_vector_int_destroy(&exend);
  splicing_vector_int_destroy(&exstart);
  SPLICING_FINALLY_CLEAN(6);
  
  return 0;
}

/* We assume that the pairs are consecutive */

int splicing_matchIso_paired(const splicing_gff_t *gff, int gene,
			     const splicing_vector_int_t *position,
			     const char **cigarstr, int readLength, 
			     int overHang,
			     const splicing_vector_t *fragmentProb,
			     int fragmentStart, double normalMean,
			     double normalVar, double numDevs,
			     splicing_matrix_t *result, 
			     splicing_matrix_int_t *fragmentLength) {
  
  size_t r, i, noreads=splicing_vector_int_size(position);
  splicing_matrix_int_t isopos;
  size_t noiso;
  size_t il;
  splicing_vector_t *myfragmentProb=(splicing_vector_t*) fragmentProb,
    vfragmentProb;

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
  
  /* Convert genomic coordinates to isoform coordinates */
  SPLICING_CHECK(splicing_matrix_int_init(&isopos, 0, 0));
  SPLICING_FINALLY(splicing_matrix_int_destroy, &isopos);
  SPLICING_CHECK(splicing_genomic_to_iso(gff, gene, position,
					 /*converter=*/ 0, &isopos));
  
  SPLICING_CHECK(splicing_matchIso(gff, gene, position, cigarstr, overHang,
				   readLength, result));
  
  if (fragmentLength) {
    SPLICING_CHECK(splicing_matrix_int_resize(fragmentLength, noiso, 
					      noreads/2));
  }
  
  for (r=0; r<noreads/2; r++) {
    for (i=0; i<noiso; i++) {
      if (MATRIX(*result, i, 2*r) && MATRIX(*result, i, 2*r+1)) {
	int frag=MATRIX(isopos, i, 2*r+1) - MATRIX(isopos, i, 2*r) +
	  readLength;
	if (frag < fragmentStart || frag >= il + fragmentStart) { 
	  MATRIX(*result, i, r) = 0.0;
	  if (fragmentLength) { MATRIX(*fragmentLength, i, r) = -1; }
	} else {
	  MATRIX(*result, i, r) = 
	    VECTOR(*myfragmentProb)[frag-fragmentStart];
	  if (fragmentLength) { MATRIX(*fragmentLength, i, r) = frag; }
	}
      } else {
	MATRIX(*result, i, r) = 0.0;
	if (fragmentLength) { MATRIX(*fragmentLength, i, r) = -1; }
      }
    }
  }    

  splicing_matrix_resize(result, noiso, noreads/2);
  
  splicing_matrix_int_destroy(&isopos);
  SPLICING_FINALLY_CLEAN(1);
  
  if (!fragmentProb) {
    splicing_vector_destroy(myfragmentProb); 
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}

int splicing_parse_cigar(const char **cigar, size_t noreads,
			 splicing_vector_int_t *numcigar,
			 splicing_vector_int_t *cigaridx, 
			 splicing_vector_int_t *cigarlength, 
			 int maxReadLength) {
  
  size_t i, pos=0;
  
  splicing_vector_int_clear(numcigar);
  SPLICING_CHECK(splicing_vector_int_resize(cigaridx, noreads+1));
  SPLICING_CHECK(splicing_vector_int_resize(cigarlength, noreads));

  for (i=0; i<noreads; i++) {
    char *s= (char*) cigar[i];
    int mode=0;			/* 0: begin, 1:middle, 2:end */
    int len=0;
    VECTOR(*cigaridx)[i] = pos;
    while (*s) {
      long l = strtol(s, &s, 10L);

      if (mode==0 && *s!='S' && *s!='H') {
	mode=1; 
      } else if (mode==1 && (*s=='S' || *s=='H')) { 
	mode=2;
      } else if (mode==2 && *s != 'S' && *s != 'H') {
	SPLICING_ERROR("Bad CIGAR string: `S' and 'H' may appear only at "
		       "the beginning and the end", SPLICING_EINVAL);
      }

      if (*s == 'M' || *s == '=') { /* MATCHING */
	if (maxReadLength > 0 && len + l > maxReadLength) { 
	  l = maxReadLength - len;
	}
	SPLICING_CHECK(splicing_vector_int_push_back(numcigar, l));
	len += l;
	pos++;
	s++;
      } else if (*s == 'N') {	/* SKIPPING */
	SPLICING_CHECK(splicing_vector_int_push_back(numcigar, -l));
	pos++;
	s++;
      } else if (*s == 'X') {	/* SEQ MISMATCH */
	if (l > 4) { SPLICING_WARNING("Long non-matching alignment"); }
	/* We count this as matching */
	if (maxReadLength > 0 && len + l > maxReadLength) { 
	  l = maxReadLength - len;
	}
	SPLICING_CHECK(splicing_vector_int_push_back(numcigar, l));
	len += l;
	pos++;
	s++;
      } else if (*s == 'S' || *s == 'H') { /* SOFT/HARD CLIPPING */
	/* We consider these 'matching' */
	if (maxReadLength > 0 && len + l > maxReadLength) { 
	  l = maxReadLength - len;
	}
	SPLICING_CHECK(splicing_vector_int_push_back(numcigar, l));
	len += l;
	pos++;
	s++;
      } else if (*s == 'D') {	/* DELETION FROM THE REFERENCE */
	if (l > 4) { SPLICING_WARNING("Long deleted alignment"); }
	/* We count this as matching */
	if (maxReadLength > 0 && len + l > maxReadLength) { 
	  l = maxReadLength - len;
	}
	SPLICING_CHECK(splicing_vector_int_push_back(numcigar, l));
	len += l;
	pos++;
	s++;
      } else if (*s == 'I') {	/* INSERTION TO THE REFERENCE */
	if (l > 4) { SPLICING_WARNING("Long inserted alignment"); }
	/* We do nothing, just ignore the part that does not appear in 
	   the genome */
	s++;
      } else {
	SPLICING_ERROR("Unsupported CIGAR string (`MNSHDI=X' are supported)", 
		       SPLICING_EINVAL);
      }
      if (maxReadLength && len == maxReadLength) { continue; }
    }
    VECTOR(*cigarlength)[i] = len;
  }
  VECTOR(*cigaridx)[i] = pos;

  return 0;
}

int splicing_solve_gene(const splicing_gff_t *gff, size_t gene, 
			int readLength, int overHang,
			const splicing_vector_int_t *position, 
			const char **cigarstr,
			splicing_matrix_t *match_matrix,
			splicing_matrix_t *assignment_matrix, 
			splicing_vector_t *expression) {

  splicing_matrix_t A;
  splicing_vector_t match;
  splicing_vector_long_t index;
  size_t no_classes;
  size_t no_reads=splicing_vector_int_size(position);
  size_t noiso;
  size_t r;
  long int nsetp;
  double rnorm;
  splicing_matrix_t *mymatch_matrix=match_matrix, vmatch_matrix;
  splicing_matrix_t *myass_matrix=assignment_matrix, vass_matrix;

  if (!match_matrix) {
    mymatch_matrix=&vmatch_matrix;
    SPLICING_CHECK(splicing_matrix_init(mymatch_matrix, 0, 0));
    SPLICING_FINALLY(splicing_matrix_destroy, mymatch_matrix);
  }

  if (!assignment_matrix) {
    myass_matrix=&vass_matrix;
    SPLICING_CHECK(splicing_matrix_init(myass_matrix, 0, 0));
    SPLICING_FINALLY(splicing_matrix_destroy, myass_matrix);
  }

  SPLICING_CHECK(splicing_assignment_matrix(gff, gene, readLength, overHang,
					    myass_matrix));
  no_classes=splicing_matrix_ncol(myass_matrix);
  noiso=splicing_matrix_nrow(myass_matrix);

  /* Calculate match vector from match matrix */
  SPLICING_CHECK(splicing_matchIso(gff, gene, position, cigarstr, overHang,
				   readLength, mymatch_matrix));
  SPLICING_CHECK(splicing_vector_init(&match, no_classes));
  SPLICING_FINALLY(splicing_vector_destroy, &match);
  for (r=0; r<no_reads; r++) {
    size_t cl, found;
    for (cl=0, found=0; !found && cl < no_classes; cl++) {
      size_t i;
      for (i=0, found=1; i<noiso && found; i++) {
	double m1=MATRIX(*mymatch_matrix, i, r);
	double m2=MATRIX(*myass_matrix, i, cl);
	found = (m1 > 0 && m2 > 0) || (m1 == 0 && m2 == 0);
      }
    }
    VECTOR(match)[cl-1] += 1;
  }

  /* This is modified, need a copy */
  SPLICING_CHECK(splicing_matrix_copy(&A, myass_matrix));
  SPLICING_FINALLY(splicing_matrix_destroy, &A);
  SPLICING_CHECK(splicing_matrix_transpose(&A));
  SPLICING_CHECK(splicing_vector_long_init(&index, 0));
  SPLICING_FINALLY(splicing_vector_long_destroy, &index);

  if (!assignment_matrix) {
    splicing_matrix_destroy(myass_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  if (!match_matrix) {
    splicing_matrix_destroy(mymatch_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  SPLICING_CHECK(splicing_nnls(&A, &match, expression, &rnorm, &index, 
			       &nsetp));
  
  splicing_vector_long_destroy(&index);
  splicing_matrix_destroy(&A);
  splicing_vector_destroy(&match);
  SPLICING_FINALLY_CLEAN(3);

  splicing_vector_scale(expression, 1.0/splicing_vector_sum(expression));

  return 0;
}

int splicing_solve_gene_paired(const splicing_gff_t *gff, size_t gene,
			       int readLength, int overHang,
			       const splicing_vector_int_t *position,
			       const char **cigarstr,
			       const splicing_vector_t *fragmentProb,
			       int fragmentStart, double normalMean,
			       double normalVar, double numDevs,
			       splicing_matrix_t *match_matrix,
			       splicing_matrix_t *assignment_matrix,
			       splicing_vector_t *expression) {
  
  splicing_matrix_t A;
  splicing_vector_t match;
  splicing_vector_long_t index;
  size_t no_classes;
  size_t no_reads=splicing_vector_int_size(position);
  size_t noiso;
  size_t r;
  long int nsetp;
  double rnorm;

  splicing_matrix_t *mymatch_matrix=match_matrix, vmatch_matrix;
  splicing_matrix_t *myass_matrix=assignment_matrix, vass_matrix;

  if (!match_matrix) {
    mymatch_matrix=&vmatch_matrix;
    SPLICING_CHECK(splicing_matrix_init(mymatch_matrix, 0, 0));
    SPLICING_FINALLY(splicing_matrix_destroy, mymatch_matrix);
  }

  if (!assignment_matrix) {
    myass_matrix=&vass_matrix;
    SPLICING_CHECK(splicing_matrix_init(myass_matrix, 0, 0));
    SPLICING_FINALLY(splicing_matrix_destroy, myass_matrix);
  }

  SPLICING_CHECK(splicing_paired_assignment_matrix(gff, gene, readLength, 
						   overHang, fragmentProb,
						   fragmentStart,
						   normalMean, normalVar,
						   numDevs,
						   myass_matrix));
  no_classes=splicing_matrix_ncol(myass_matrix);
  noiso=splicing_matrix_nrow(myass_matrix);

  /* Calculate match vector from match matrix */
  SPLICING_CHECK(splicing_matchIso_paired(gff, gene, position, cigarstr, 
					  readLength, overHang, fragmentProb,
					  fragmentStart, normalMean, 
					  normalVar, numDevs, mymatch_matrix, 
					  0));
  SPLICING_CHECK(splicing_vector_init(&match, no_classes));
  SPLICING_FINALLY(splicing_vector_destroy, &match);
  for (r=0; r<no_reads/2; r++) {
    size_t cl, found;
    for (cl=0, found=0; !found && cl < no_classes; cl++) {
      size_t i;
      for (i=0, found=1; i<noiso && found; i++) {
	double m1=MATRIX(*mymatch_matrix, i, r);
	double m2=MATRIX(*myass_matrix, i, cl);
	found = (m1 > 0 && m2 > 0) || (m1 == 0 && m2 == 0);
      }
    }
    if (!found) {
      SPLICING_WARNING("Read does not match any assignment class");
      /* SPLICING_ERROR("Read does not match any assignment class",  */
      /* 		     SPLICING_EINTERNAL); */
    }
    VECTOR(match)[cl-1] += 1;
  }

  /* This is modified, need a copy */
  SPLICING_CHECK(splicing_matrix_copy(&A, myass_matrix));
  SPLICING_FINALLY(splicing_matrix_destroy, &A);
  SPLICING_CHECK(splicing_matrix_transpose(&A));
  SPLICING_CHECK(splicing_vector_long_init(&index, 0));
  SPLICING_FINALLY(splicing_vector_long_destroy, &index);

  if (!assignment_matrix) {
    splicing_matrix_destroy(myass_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  if (!match_matrix) {
    splicing_matrix_destroy(mymatch_matrix);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  SPLICING_CHECK(splicing_nnls(&A, &match, expression, &rnorm, &index,
			       &nsetp));
  
  splicing_vector_long_destroy(&index);
  splicing_matrix_destroy(&A);
  splicing_vector_destroy(&match);
  SPLICING_FINALLY_CLEAN(3);

  splicing_vector_scale(expression, 1.0/splicing_vector_sum(expression));

  return 0;
  
}
