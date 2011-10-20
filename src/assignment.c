
#include "splicing.h"
#include "splicing_error.h"

#include <stdlib.h>
#include <string.h>

typedef struct {
  splicing_vector_int_t *mp;
  splicing_vector_int_t *mppos;
} splicing_i_assignmat_data_t;

int splicing_i_assignmat_cmp(void *pdata, const void *a, const void *b) {
  const splicing_i_assignmat_data_t *data =
    (splicing_i_assignmat_data_t*) pdata;
  const splicing_vector_int_t *mp=data->mp;
  const splicing_vector_int_t *mppos=data->mppos;
  int *aa=(int*) a;
  int *bb=(int*) b;
  int p1=VECTOR(*mppos)[*aa]-1;
  int p2=VECTOR(*mppos)[*bb]-1;
  
  if (p1 < 0 && p2 < 0) { return 0; }
  if (p1 < 0) { return -1; }
  if (p2 < 0) { return 1; }
  
  while (VECTOR(*mp)[p1] != 0 && VECTOR(*mp)[p2] != 0 && 
	 VECTOR(*mp)[p1] == VECTOR(*mp)[p2]) { p1++; p2++; }
  
  if (VECTOR(*mp)[p1] == 0 && VECTOR(*mp)[p2] == 0) { return 0; }
  if (VECTOR(*mp)[p1] == 0) { return -1; } 
  if (VECTOR(*mp)[p2] == 0) { return  1; }
  return VECTOR(*mp)[p1] < VECTOR(*mp)[p2] ? -1 : 1;
}

int splicing_i_assignmat_cmp2(void *data, const void *a, const void *b) {
  splicing_matrix_t *m=(splicing_matrix_t*) data;
  int aa=*(int*)a, bb=*(int*)b;
  double *aaa=&MATRIX(*m, 0, aa);
  double *bbb=&MATRIX(*m, 0, bb);
  int i, nrow=splicing_matrix_nrow(m);
  
  for (i=0; i<nrow; i++) {
    if (aaa[i]==0 && bbb[i] != 0) { return -1; }
    if (aaa[i]!=0 && bbb[i] == 0) { return  1; }
  }
  return 0;
}

int splicing_i_assignmat_simplify(splicing_matrix_t *mat) {
  splicing_matrix_t tmp;
  splicing_vector_int_t order;
  int i, j;
  int ncol=splicing_matrix_ncol(mat);
  int nrow=splicing_matrix_nrow(mat);
  
  SPLICING_CHECK(splicing_matrix_copy(&tmp, mat));
  SPLICING_FINALLY(splicing_matrix_destroy, &tmp);
  SPLICING_CHECK(splicing_vector_int_init(&order, ncol));
  SPLICING_FINALLY(splicing_vector_int_destroy, &order);
  for (i=0; i<ncol; i++) { VECTOR(order)[i]=i; }
  splicing_qsort_r(VECTOR(order), ncol, sizeof(int), (void*) &tmp, 
		   splicing_i_assignmat_cmp2);

  memcpy(&MATRIX(*mat, 0, 0), &MATRIX(tmp, 0, VECTOR(order)[0]), 
	 sizeof(double) * nrow);
  for (i=1, j=0; i<ncol; i++) {
    if (!splicing_i_assignmat_cmp2(&tmp, &(VECTOR(order)[i-1]), 
				   &(VECTOR(order)[i]))) {
      int k;
      for (k=0; k<nrow; k++) {
	MATRIX(*mat, k, j) += MATRIX(tmp, k, VECTOR(order)[i]);
      }
    } else {
      j++;
      memcpy(&MATRIX(*mat, 0, j), &MATRIX(tmp, 0, VECTOR(order)[i]), 
	     sizeof(double) * nrow);
    }
  }
  SPLICING_CHECK(splicing_matrix_resize(mat, nrow, j+1));

  splicing_vector_int_destroy(&order);
  splicing_matrix_destroy(&tmp);
  SPLICING_FINALLY_CLEAN(2);

  return 0;
}

int splicing_assignment_matrix(const splicing_gff_t *gff, size_t gene,
			       int readLength, int overHang, 
			       splicing_matrix_t *matrix) {
  size_t noiso;
  splicing_vector_int_t exstart, exend, exidx;
  size_t genestart=VECTOR(gff->start)[ VECTOR(gff->genes)[gene] ];
  size_t geneend=VECTOR(gff->end)[ VECTOR(gff->genes)[gene] ];
  size_t i, elen;
  size_t p=0, lastp=geneend - genestart - readLength + 1;
  splicing_vector_int_t cigar, mp, mppos;
  splicing_vector_int_t isoseq, isomatch;
  splicing_i_assignmat_data_t mpdata = { &mp, &mppos };

  if (overHang > 1) { 
    SPLICING_ERROR("Overhang is not implemented in assignment matrix yet.",
		   SPLICING_UNIMPLEMENTED);
  }

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));

  SPLICING_CHECK(splicing_vector_int_init(&cigar, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &cigar);
  SPLICING_CHECK(splicing_vector_int_init(&exstart, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exstart);
  SPLICING_CHECK(splicing_vector_int_init(&exend, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exend);
  SPLICING_CHECK(splicing_vector_int_init(&exidx, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exidx);
  SPLICING_CHECK(splicing_gff_exon_start_end(gff, &exstart, &exend, &exidx,
					     gene));
  elen=splicing_vector_int_size(&exstart);
  for (i=0; i<elen; i++) {
    VECTOR(exstart)[i] -= genestart;
    VECTOR(exend)[i] -= genestart;
  }
  
  SPLICING_CHECK(splicing_numeric_cigar(&exstart, &exend, &exidx, noiso, 0,
					&cigar, /* start= */ 0));

  splicing_vector_int_destroy(&exidx);
  splicing_vector_int_destroy(&exend);
  splicing_vector_int_destroy(&exstart);
  SPLICING_FINALLY_CLEAN(3);
  
  SPLICING_CHECK(splicing_vector_int_init(&mp, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &mp);
  SPLICING_CHECK(splicing_vector_int_init(&mppos, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &mppos);
  SPLICING_CHECK(splicing_vector_int_init(&isoseq, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isoseq);
  SPLICING_CHECK(splicing_vector_int_init(&isomatch, noiso));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isomatch);

  SPLICING_CHECK(splicing_matrix_resize(matrix, noiso, 0));

  while (p <= lastp) {
    size_t pos=0, nextp, tmppos, ipos, l;

    /* Calculate the initial CIGAR strings for all isoforms */
    splicing_vector_int_clear(&mp);
    SPLICING_CHECK(splicing_vector_int_push_back(&mp, 0)); l=1;
    for (i=0; i<noiso; i++) {
      int rl=readLength;
      pos++;
      VECTOR(mppos)[i]=l+1;
      if (VECTOR(cigar)[pos] <= 0) { 
	VECTOR(mppos)[i]=0;
	while (VECTOR(cigar)[pos] != 0) { pos++; }
      } else {
	while (VECTOR(cigar)[pos] != 0) {
	  if (VECTOR(cigar)[pos] < 0) {
	    SPLICING_CHECK(splicing_vector_int_push_back(&mp, 
					 VECTOR(cigar)[pos])); l++;
	    pos++;
	  } else if (rl <= VECTOR(cigar)[pos]) {
	    SPLICING_CHECK(splicing_vector_int_push_back(&mp, rl)); l++;
	    while (VECTOR(cigar)[pos] != 0) { pos++; }
	    rl = 0;
	    break;
	  } else {
	    SPLICING_CHECK(splicing_vector_int_push_back(&mp, 
					 VECTOR(cigar)[pos])); l++;
	    rl=rl-VECTOR(cigar)[pos];
	    pos++;
	  }
	}
	if (rl != 0) { VECTOR(mppos)[i] = 0; }
      }
      SPLICING_CHECK(splicing_vector_int_push_back(&mp, 0)); l++;
    } /* i<noiso */

    /* Calculate the next position to check */
    nextp = lastp + 1 - p;
    tmppos = 0;
    for (i=0; i<noiso; i++) {
      tmppos++;
      ipos=tmppos;
      if (VECTOR(cigar)[tmppos]==0) { continue; }
      if (abs(VECTOR(cigar)[tmppos]) < nextp) { 
	nextp = abs(VECTOR(cigar)[tmppos]);
      }
      if (VECTOR(cigar)[tmppos] > 0) {
	size_t pnextp=nextp;
	int j=ipos;
	int rl2=readLength;
	while (VECTOR(cigar)[j] != 0) {
	  if (VECTOR(cigar)[j] >= rl2) { 
	    pnextp=VECTOR(cigar)[j]-rl2+1; 
	    break;
	  } else if (VECTOR(cigar)[j] > 0) {
	    rl2 = rl2-VECTOR(cigar)[j];
	  }
	  j++;
	}
	if (pnextp < nextp) { nextp = pnextp; }
      }
      while (VECTOR(cigar)[tmppos] != 0) { tmppos++; }
    }

    /* Record the assignment classes. 
       First we sort the cigar prefixes; then we calculate the
       assignment class of the current position. Finally, we check
       whether this class is already in the result and then update/add
       it. */
    for (i=0; i<noiso; i++) { VECTOR(isoseq)[i]=i; }
    splicing_qsort_r(VECTOR(isoseq), noiso, sizeof(int), &mpdata, 
		     splicing_i_assignmat_cmp);
    
    for (i=0; i<noiso && VECTOR(mppos)[ VECTOR(isoseq)[i] ]==0; i++) ;
    while (i<noiso) {
      int col, j;
      i++;
      splicing_vector_int_null(&isomatch);
      VECTOR(isomatch)[ VECTOR(isoseq)[i-1] ] = 1;
      while (i < noiso &&  ! splicing_i_assignmat_cmp(&mpdata, 
					      &(VECTOR(isoseq)[i-1]),
					      &(VECTOR(isoseq)[i]))) {
	VECTOR(isomatch)[ VECTOR(isoseq)[i] ] = 1;
	i++;
      }
      SPLICING_CHECK(splicing_matrix_add_cols(matrix, 1));
      col=splicing_matrix_ncol(matrix)-1;
      for (j=0; j<noiso; j++) { 
	MATRIX(*matrix, j, col) = VECTOR(isomatch)[j] * nextp;
      }
    }

    /* Update the CIGAR strings */
    tmppos=0, ipos=0;
    for (i=0; i<noiso; i++) {
      VECTOR(cigar)[ipos] = 0;
      tmppos++; ipos++;
      if (VECTOR(cigar)[tmppos] != 0) {
	VECTOR(cigar)[ipos] = 
	  (VECTOR(cigar)[tmppos] < 0 ? -1 : 1) * 
	  (abs(VECTOR(cigar)[tmppos])-nextp);
	tmppos++;
	if (VECTOR(cigar)[ipos] !=0 ) ipos++;
      }
      while (VECTOR(cigar)[tmppos] != 0) {
	VECTOR(cigar)[ipos]=VECTOR(cigar)[tmppos];
	tmppos++;
	ipos++;
      }      
    }
    VECTOR(cigar)[ipos] = 0;
    
    p = p + nextp;


  } /* p <= lastp */
  
  splicing_vector_int_destroy(&isomatch);
  splicing_vector_int_destroy(&isoseq);
  splicing_vector_int_destroy(&mppos);
  splicing_vector_int_destroy(&mp);
  splicing_vector_int_destroy(&cigar);
  SPLICING_FINALLY_CLEAN(5);

  splicing_i_assignmat_simplify(matrix);
  
  return 0;
}

/* TODO: speedup, many things are recalculated in the 
   splicing_assignment_matrix calls. */

int splicing_paired_assignment_matrix(const splicing_gff_t *gff, size_t gene,
				      int readLength, int overHang,
				      const splicing_vector_t *fragmentProb,
				      int fragmentStart, double normalMean,
				      double normalVar, double numDevs,
				      splicing_matrix_t *matrix) {
  size_t noiso;
  int i, il;
  splicing_matrix_t tmpmat;
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

  SPLICING_CHECK(splicing_matrix_init(&tmpmat, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &tmpmat);
  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));
  
  SPLICING_CHECK(splicing_matrix_resize(matrix, noiso, 0));
  for (i=0; i<il; i++) {

    int c, tmpncol, m, mncol=splicing_matrix_ncol(matrix);
    int myrl=i + fragmentStart;
    double fact=VECTOR(*myfragmentProb)[i];

    SPLICING_CHECK(splicing_assignment_matrix(gff, gene, myrl, overHang, 
					      &tmpmat));
    tmpncol=splicing_matrix_ncol(&tmpmat);
    for (c=0; c<tmpncol; c++) {

      int j, found;
      /* find the corresponding column in 'matrix', if it exists */
      for (m=0, found=0; !found && m<mncol; m++) {
	for (j=0, found=1; found && j<noiso; j++) {
	  double m1=MATRIX(tmpmat, j, c);
	  double m2=MATRIX(*matrix, j, m);
	  found = (m1 > 0 && m2 > 0) || (m1==0 && m2==0);
	}
      }

      if (found) {
	for (j=0; j<noiso; j++) {
	  MATRIX(*matrix, j, m-1) += fact * MATRIX(tmpmat, j, c);
	}
      } else {
	SPLICING_CHECK(splicing_matrix_add_cols(matrix, 1));
	m=splicing_matrix_ncol(matrix)-1;
	for (j=0; j<noiso; j++) { 
	  MATRIX(*matrix, j, m) = fact * MATRIX(tmpmat, j, c);
	}
      }

    } /* c < tmpncol */
  } /* i < il */

  splicing_matrix_destroy(&tmpmat);
  SPLICING_FINALLY_CLEAN(1);

  if (!fragmentProb) { 
    splicing_vector_destroy(myfragmentProb); 
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}
