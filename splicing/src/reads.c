
#include <ctype.h>

#include "splicing.h"
#include "splicing_error.h"

#include "sam.h"

void splicing_i_bam_destroy1(bam1_t *read) {
  bam_destroy1(read);
}

typedef struct splicing_i_read_sambam_sort_data_t {
  splicing_vector_int_t *chr;
  splicing_vector_int_t *pos;
  splicing_vector_int_t *flag;
  splicing_vector_int_t *pairpos;
} splicing_i_read_sambam_sort_data_t;

int splicing_i_cmp_reads(void *pdata, const void *a, const void *b) {
  const splicing_i_read_sambam_sort_data_t *data = 
    (splicing_i_read_sambam_sort_data_t *) pdata;
  int aa = *(int*) a, bb = *(int*) b;
  int apos, bpos, arpos, brpos;

  if (VECTOR(*data->chr)[aa] < VECTOR(*data->chr)[bb]) { 
    return -1; 
  } else if (VECTOR(*data->chr)[aa] > VECTOR(*data->chr)[bb]) { 
    return 1;
  }
  
  apos=VECTOR(*data->pos)[aa];
  bpos=VECTOR(*data->pos)[bb];

  arpos=VECTOR(*data->pairpos)[aa]; if (arpos > apos) { arpos=apos; }
  brpos=VECTOR(*data->pairpos)[bb]; if (brpos > bpos) { brpos=bpos; }

  if (arpos < brpos) { return -1; } else if (arpos > brpos) { return 1; }
  if (apos  < bpos ) { return -1; } else if (apos  > bpos ) { return 1; }
  
  return 0;
}

/* TODO: use qname as well for pairing */

int splicing_read_sambam(const char *filename, 
			 splicing_strvector_t *chrname, 
			 splicing_vector_int_t *chrlen,
			 splicing_vector_int_t *chr,
			 splicing_strvector_t *qname,
			 splicing_strvector_t *cigar,
			 splicing_vector_int_t *position,
			 splicing_vector_int_t *flag,
			 splicing_vector_int_t *pairpos,
			 int *noPairs, int *noSingles, int *paired,
			 splicing_vector_int_t *mapq,
			 splicing_vector_int_t *rnext,
			 splicing_vector_int_t *tlen,
			 splicing_strvector_t *seq,
			 splicing_strvector_t *qual) {

  samfile_t *infile=0;
  int bytesread;
  bam1_t *read = bam_init1();
  char buffer[4096];
  char *cigarcode="MIDNSHP";
  int i;

  SPLICING_FINALLY(splicing_i_bam_destroy1, read);
  
  splicing_strvector_clear(chrname);
  splicing_vector_int_clear(chrlen);
  splicing_vector_int_clear(chr);
  splicing_strvector_clear(qname);
  splicing_strvector_clear(cigar);
  splicing_vector_int_clear(position);
  splicing_vector_int_clear(flag);
  splicing_vector_int_clear(pairpos);
  splicing_vector_int_clear(mapq);
  splicing_vector_int_clear(rnext);
  splicing_vector_int_clear(tlen);
  splicing_strvector_clear(seq);
  splicing_strvector_clear(qual);
  *noPairs = *noSingles = 0;

  infile=samopen(filename, "r", /*aux=*/ 0);
  if (!infile) { 
    SPLICING_ERROR("Cannot open SAM/BAM file", SPLICING_EFILE);
  }

  SPLICING_CHECK(splicing_vector_int_resize(chrlen,
					    infile->header->n_targets));
  for (i=0; i<infile->header->n_targets; i++) {
    char *s=infile->header->target_name[i];
    SPLICING_CHECK(splicing_strvector_append(chrname, s));
    VECTOR(*chrlen)[i] = infile->header->target_len[i];
  }

  while ( (bytesread=samread(infile, read)) >= 0) {
    int i, ncigar=read->core.n_cigar;
    char *bufptr=buffer;
    uint32_t *actcigar=bam1_cigar(read);
    uint8_t *s = bam1_seq(read), *t = bam1_qual(read);

    SPLICING_CHECK(splicing_vector_int_push_back(position,
						 read->core.pos+1));
    SPLICING_CHECK(splicing_strvector_append2(qname, 
					      bam1_qname(read), 
					      read->core.l_qname-1));
    if (read->core.n_cigar == 0) {
      SPLICING_CHECK(splicing_strvector_append(cigar, "*"));
    } else {
      for (i=0; i<ncigar; i++) {
	int l=snprintf(bufptr, sizeof(buffer)-(bufptr-buffer)*sizeof(char),
		       "%i%c", (int) (actcigar[i] >> BAM_CIGAR_SHIFT), 
		       cigarcode[actcigar[i] & BAM_CIGAR_MASK]);
	bufptr += l;
      }
      SPLICING_CHECK(splicing_strvector_append2(cigar, buffer,
						bufptr-buffer));
    }
    
    SPLICING_CHECK(splicing_vector_int_push_back(chr, read->core.tid));
    SPLICING_CHECK(splicing_vector_int_push_back(pairpos, 
						 read->core.mpos+1));
    SPLICING_CHECK(splicing_vector_int_push_back(flag, read->core.flag));
    SPLICING_CHECK(splicing_vector_int_push_back(rnext, read->core.mtid));
    SPLICING_CHECK(splicing_vector_int_push_back(mapq, read->core.qual));
    SPLICING_CHECK(splicing_vector_int_push_back(tlen, read->core.isize));
    
    if (!read->core.l_qseq) {
      SPLICING_CHECK(splicing_strvector_append(seq, "*"));
    } else {
      bufptr=buffer;
      for (i=0; i<read->core.l_qseq; i++) { 
	*bufptr = bam_nt16_rev_table[bam1_seqi(s, i)];
	bufptr++;
      }
      SPLICING_CHECK(splicing_strvector_append2(seq, buffer, bufptr-buffer));
    }

    if (!read->core.l_qseq || t[0] == 0xff) {
      SPLICING_CHECK(splicing_strvector_append(qual, "*"));
    } else {
      bufptr=buffer;
      for (i=0; i<read->core.l_qseq; i++) { 
	*bufptr = t[i] + 33;
	bufptr++;
      }
      SPLICING_CHECK(splicing_strvector_append2(qual, buffer, bufptr-buffer));
    }

    if (read->core.mpos < 0) {
      *noSingles  += 1;
    } else {
      *noPairs += 1;
    }
  }

  *noPairs /= 2;
  *paired = *noPairs ? 1 : 0;

  if (bytesread < -1) {
    SPLICING_WARNING("Truncated SAM/BAM file");
  }

  /* We need to find the pairs of the paired-end reads, if there is any.
     For this we first order the reads according to their RNAME
     (chromosome), plus their position. Then we search for the second
     pair of each first pair. */

  if (*noPairs != 0) {
    splicing_vector_char_t taken;
    splicing_vector_int_t idx;
    splicing_vector_t idx2;
    size_t i, pos, noReads=(*noPairs) * 2 + (*noSingles);
    splicing_i_read_sambam_sort_data_t data = 
      { chr, position, flag, pairpos };

    SPLICING_CHECK(splicing_vector_int_init(&idx, noReads));
    SPLICING_FINALLY(splicing_vector_int_destroy, &idx);
    SPLICING_CHECK(splicing_vector_init(&idx2, noReads));
    SPLICING_FINALLY(splicing_vector_destroy, &idx2);
    SPLICING_CHECK(splicing_vector_char_init(&taken, noReads));
    SPLICING_FINALLY(splicing_vector_char_destroy, &taken);
    for (i=0; i<noReads; i++) { VECTOR(idx)[i] = i; }
    splicing_qsort_r(VECTOR(idx), noReads, sizeof(int), (void*) &data, 
		     splicing_i_cmp_reads);
    
    for (i=0, pos=0; i<noReads; i++) {
      int curr=VECTOR(idx)[i];
      if (VECTOR(*pairpos)[curr] == 0) { 
	/* has no pair, add it */
	VECTOR(idx2)[pos++] = curr;
      } else if (VECTOR(*pairpos)[curr] < VECTOR(*position)[curr])  {
	/* read and pair already counted, do nothing */
      } else {
	/* search for pair */
	int needle=VECTOR(*pairpos)[curr], ppos, found;
	for (ppos=i+1, found=0; ! found && ppos < noReads; ppos++) {
	  int rp=VECTOR(idx)[ppos];
	  found = VECTOR(*position)[rp] == needle && 
	    VECTOR(*pairpos)[rp] == VECTOR(*position)[curr] &&
	    0x80 & VECTOR(*flag)[rp] && 
	    ! VECTOR(taken)[rp];
	}
	if (!found) { 
	  SPLICING_WARNING("Fragment missing from read pair");
	} else { 
	  int rp=VECTOR(idx)[ppos-1];
	  VECTOR(idx2)[pos++] = curr;
	  VECTOR(idx2)[pos++] = rp;
	  VECTOR(taken)[rp] = 1;
	}
      }
    }

    splicing_vector_char_destroy(&taken);
    SPLICING_FINALLY_CLEAN(1);

    /* We have the correct order now, reorder the vectors */
    splicing_vector_int_iindex(chr, &idx2);
    splicing_strvector_permute(qname, &idx2);
    splicing_strvector_permute(cigar, &idx2);
    splicing_vector_int_iindex(position, &idx2);
    splicing_vector_int_iindex(flag, &idx2);
    splicing_vector_int_iindex(pairpos, &idx2);

    splicing_vector_destroy(&idx2);
    splicing_vector_int_destroy(&idx);
    SPLICING_FINALLY_CLEAN(2);
  }

  bam_destroy1(read);
  SPLICING_FINALLY_CLEAN(1);

  samclose(infile);
  
  return 0;
}
