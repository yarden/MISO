
#include <ctype.h>

#include "splicing.h"
#include "splicing_error.h"

#include "sam.h"

int splicing_reads_init(splicing_reads_t *reads) {
  reads->noPairs = reads->noSingles = 0;
  reads->paired = 0;

  SPLICING_CHECK(splicing_strvector_init(&reads->chrname, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->chrname);
  SPLICING_CHECK(splicing_vector_int_init(&reads->chrlen, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->chrlen);
  SPLICING_CHECK(splicing_vector_int_init(&reads->chr, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->chr);
  SPLICING_CHECK(splicing_strvector_init(&reads->qname, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->qname);
  SPLICING_CHECK(splicing_strvector_init(&reads->cigar, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->cigar);
  SPLICING_CHECK(splicing_vector_int_init(&reads->position, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->position);
  SPLICING_CHECK(splicing_vector_int_init(&reads->flags, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->flags);
  SPLICING_CHECK(splicing_vector_int_init(&reads->pairpos, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->pairpos);
  SPLICING_CHECK(splicing_vector_int_init(&reads->mapq, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->mapq);
  SPLICING_CHECK(splicing_vector_int_init(&reads->rnext, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->rnext);
  SPLICING_CHECK(splicing_vector_int_init(&reads->tlen, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->tlen);
  SPLICING_CHECK(splicing_strvector_init(&reads->seq, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->seq);
  SPLICING_CHECK(splicing_strvector_init(&reads->qual, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->qual);
  
  SPLICING_FINALLY_CLEAN(13);

  return 0;
}

void splicing_reads_destroy(splicing_reads_t *reads) {
  splicing_strvector_destroy(&reads->chrname);
  splicing_vector_int_destroy(&reads->chrlen);
  splicing_vector_int_destroy(&reads->chr);
  splicing_strvector_destroy(&reads->qname);
  splicing_strvector_destroy(&reads->cigar);
  splicing_vector_int_destroy(&reads->position);
  splicing_vector_int_destroy(&reads->flags);
  splicing_vector_int_destroy(&reads->pairpos);
  splicing_vector_int_destroy(&reads->mapq);
  splicing_vector_int_destroy(&reads->rnext);
  splicing_vector_int_destroy(&reads->tlen);
  splicing_strvector_destroy(&reads->seq);
  splicing_strvector_destroy(&reads->qual);
}

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
			 splicing_reads_t *reads) {

  samfile_t *infile=0;
  int bytesread;
  bam1_t *read = bam_init1();
  char buffer[4096];
  char *cigarcode="MIDNSHP";
  int i;

  SPLICING_FINALLY(splicing_i_bam_destroy1, read);
  
  splicing_strvector_clear(&reads->chrname);
  splicing_vector_int_clear(&reads->chrlen);
  splicing_vector_int_clear(&reads->chr);
  splicing_strvector_clear(&reads->qname);
  splicing_strvector_clear(&reads->cigar);
  splicing_vector_int_clear(&reads->position);
  splicing_vector_int_clear(&reads->flags);
  splicing_vector_int_clear(&reads->pairpos);
  splicing_vector_int_clear(&reads->mapq);
  splicing_vector_int_clear(&reads->rnext);
  splicing_vector_int_clear(&reads->tlen);
  splicing_strvector_clear(&reads->seq);
  splicing_strvector_clear(&reads->qual);
  reads->noPairs = reads->noSingles = 0;

  infile=samopen(filename, "r", /*aux=*/ 0);
  if (!infile) { 
    SPLICING_ERROR("Cannot open SAM/BAM file", SPLICING_EFILE);
  }

  SPLICING_CHECK(splicing_vector_int_resize(&reads->chrlen,
					    infile->header->n_targets));
  for (i=0; i<infile->header->n_targets; i++) {
    char *s=infile->header->target_name[i];
    SPLICING_CHECK(splicing_strvector_append(&reads->chrname, s));
    VECTOR(reads->chrlen)[i] = infile->header->target_len[i];
  }

  while ( (bytesread=samread(infile, read)) >= 0) {
    int i, ncigar=read->core.n_cigar;
    char *bufptr=buffer;
    uint32_t *actcigar=bam1_cigar(read);
    uint8_t *s = bam1_seq(read), *t = bam1_qual(read);

    SPLICING_CHECK(splicing_vector_int_push_back(&reads->position,
						 read->core.pos+1));
    SPLICING_CHECK(splicing_strvector_append2(&reads->qname, 
					      bam1_qname(read), 
					      read->core.l_qname-1));
    if (read->core.n_cigar == 0) {
      SPLICING_CHECK(splicing_strvector_append(&reads->cigar, "*"));
    } else {
      for (i=0; i<ncigar; i++) {
	int l=snprintf(bufptr, sizeof(buffer)-(bufptr-buffer)*sizeof(char),
		       "%i%c", (int) (actcigar[i] >> BAM_CIGAR_SHIFT), 
		       cigarcode[actcigar[i] & BAM_CIGAR_MASK]);
	bufptr += l;
      }
      SPLICING_CHECK(splicing_strvector_append2(&reads->cigar, buffer,
						bufptr-buffer));
    }
    
    SPLICING_CHECK(splicing_vector_int_push_back(&reads->chr,
						 read->core.tid));
    SPLICING_CHECK(splicing_vector_int_push_back(&reads->pairpos, 
						 read->core.mpos+1));
    SPLICING_CHECK(splicing_vector_int_push_back(&reads->flags, 
						 read->core.flag));
    SPLICING_CHECK(splicing_vector_int_push_back(&reads->rnext, 
						 read->core.mtid));
    SPLICING_CHECK(splicing_vector_int_push_back(&reads->mapq, 
						 read->core.qual));
    SPLICING_CHECK(splicing_vector_int_push_back(&reads->tlen,
						 read->core.isize));
    
    if (!read->core.l_qseq) {
      SPLICING_CHECK(splicing_strvector_append(&reads->seq, "*"));
    } else {
      bufptr=buffer;
      for (i=0; i<read->core.l_qseq; i++) { 
	*bufptr = bam_nt16_rev_table[bam1_seqi(s, i)];
	bufptr++;
      }
      SPLICING_CHECK(splicing_strvector_append2(&reads->seq, buffer, 
						bufptr-buffer));
    }

    if (!read->core.l_qseq || t[0] == 0xff) {
      SPLICING_CHECK(splicing_strvector_append(&reads->qual, "*"));
    } else {
      bufptr=buffer;
      for (i=0; i<read->core.l_qseq; i++) { 
	*bufptr = t[i] + 33;
	bufptr++;
      }
      SPLICING_CHECK(splicing_strvector_append2(&reads->qual, buffer,
						bufptr-buffer));
    }

    if (read->core.mpos < 0) {
      reads->noSingles  += 1;
    } else {
      reads->noPairs += 1;
    }
  }

  reads->noPairs /= 2;
  reads->paired = reads->noPairs ? 1 : 0;

  if (bytesread < -1) {
    SPLICING_WARNING("Truncated SAM/BAM file");
  }

  /* We need to find the pairs of the paired-end reads, if there is any.
     For this we first order the reads according to their RNAME
     (chromosome), plus their position. Then we search for the second
     pair of each first pair. */

  if (reads->noPairs != 0) {
    splicing_vector_char_t taken;
    splicing_vector_int_t idx;
    splicing_vector_t idx2;
    size_t i, pos, noReads=(reads->noPairs) * 2 + (reads->noSingles);
    splicing_i_read_sambam_sort_data_t data = 
      { &reads->chr, &reads->position, &reads->flags, &reads->pairpos };

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
      if (VECTOR(reads->pairpos)[curr] == 0) { 
	/* has no pair, add it */
	VECTOR(idx2)[pos++] = curr;
      } else if (VECTOR(reads->pairpos)[curr] < 
		 VECTOR(reads->position)[curr])  {
	/* read and pair already counted, do nothing */
      } else {
	/* search for pair */
	int needle=VECTOR(reads->pairpos)[curr], ppos, found;
	for (ppos=i+1, found=0; ! found && ppos < noReads; ppos++) {
	  int rp=VECTOR(idx)[ppos];
	  found = VECTOR(reads->position)[rp] == needle && 
	    VECTOR(reads->pairpos)[rp] == VECTOR(reads->position)[curr] &&
	    0x80 & VECTOR(reads->flags)[rp] && 
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
    splicing_vector_int_iindex(&reads->chr, &idx2);
    splicing_strvector_permute(&reads->qname, &idx2);
    splicing_strvector_permute(&reads->cigar, &idx2);
    splicing_vector_int_iindex(&reads->position, &idx2);
    splicing_vector_int_iindex(&reads->flags, &idx2);
    splicing_vector_int_iindex(&reads->pairpos, &idx2);

    splicing_vector_destroy(&idx2);
    splicing_vector_int_destroy(&idx);
    SPLICING_FINALLY_CLEAN(2);
  }

  bam_destroy1(read);
  SPLICING_FINALLY_CLEAN(1);

  samclose(infile);
  
  return 0;
}
