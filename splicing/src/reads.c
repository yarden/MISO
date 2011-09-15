
#include <ctype.h>
#include <string.h>

#include "splicing.h"
#include "splicing_error.h"

#include "sam.h"
#include "Rsplicing.h"

bam_index_t *bam_index_load_core(FILE *fp);

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
  SPLICING_CHECK(splicing_vector_int_init(&reads->mypair, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->mypair);
  SPLICING_CHECK(splicing_strvector_init(&reads->attributes, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->attributes);
  
  SPLICING_FINALLY_CLEAN(15);

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
  splicing_vector_int_destroy(&reads->mypair);
  splicing_strvector_destroy(&reads->attributes);
}

int splicing_reads_clear(splicing_reads_t *reads) {
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
  splicing_vector_int_clear(&reads->mypair);
  splicing_strvector_clear(&reads->attributes);
  reads->noPairs = reads->noSingles = 0;
  return 0;
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
  int apos, bpos, arpos, brpos, aapos, bbpos;

  if (VECTOR(*data->chr)[aa] < VECTOR(*data->chr)[bb]) { 
    return -1; 
  } else if (VECTOR(*data->chr)[aa] > VECTOR(*data->chr)[bb]) { 
    return 1;
  }
  
  aapos=apos=VECTOR(*data->pos)[aa];
  bbpos=bpos=VECTOR(*data->pos)[bb];

  arpos=VECTOR(*data->pairpos)[aa]; 
  if (arpos > apos) { int tmp=arpos; arpos=apos; apos=tmp; } 
  brpos=VECTOR(*data->pairpos)[bb];
  if (brpos > bpos) { int tmp=brpos; brpos=bpos; bpos=tmp; } 

  if (arpos < brpos) { return -1; } else if (arpos > brpos) { return 1; }
  if (apos  < bpos ) { return -1; } else if (apos  > bpos ) { return 1; }
  if (aapos < bbpos) { return -1; } else if (aapos > bbpos) { return 1; }  
  
  return 0;
}

#define REMAINING (sizeof(buffer)-(bufptr-buffer)*sizeof(char))

int splicing_i_add_read(splicing_reads_t *reads, const bam1_t *read) {
  int i, ncigar=read->core.n_cigar;
  char buffer[4096];
  char *cigarcode="MIDNSHP";
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
      int l=snprintf(bufptr, REMAINING,
		     "%i%c", (int) (actcigar[i] >> BAM_CIGAR_SHIFT), 
		     cigarcode[actcigar[i] & BAM_CIGAR_MASK]);
      bufptr += l;
    }
    SPLICING_CHECK(splicing_strvector_append2(&reads->cigar, buffer,
					      bufptr-buffer));
  }
  
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->chr, read->core.tid));
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

  SPLICING_CHECK(splicing_vector_int_push_back(&reads->mypair, -1));

  /* Extra columns */
  
  bufptr = buffer;
  s = bam1_aux(read);
  i=0;
  while (s < read->data + read->data_len) {
    uint8_t type, key[2];
    int l;
    key[0] = s[0]; key[1] = s[1];
    s += 2; type = *s; ++s;
    if (i!=0) { *bufptr = '\t'; bufptr++; }
    l=snprintf(bufptr, REMAINING, "%c%c:", (int) key[0], (int) key[1]);
    bufptr += l;
    switch (type) { 
      uint8_t sub_type;
      int32_t n;
      int j;
    case 'A': 
      l=snprintf(bufptr, REMAINING, "A:%c", (int) *s); ++s; bufptr += l;
      break;
    case 'C':
      l=snprintf(bufptr, REMAINING, "i:%u", *(uint8_t*)s); ++s; bufptr += l;
      break;
    case 'c':
      l=snprintf(bufptr, REMAINING, "i:%i", *(int8_t*)s); ++s; bufptr += l;
      break;
    case 'S':
      l=snprintf(bufptr, REMAINING, "i:%u", *(uint16_t*)s); s+=2; bufptr += l;
      break;
    case 's': 
      l=snprintf(bufptr, REMAINING, "i:%i", *(int16_t*)s); s+=2; bufptr += l;
      break;
    case 'I':
      l=snprintf(bufptr, REMAINING, "i:%u", *(uint32_t*)s); s+=4; bufptr +=l;
      break;
    case 'i':
      l=snprintf(bufptr, REMAINING, "i:%i", *(int32_t*)s); s+=4; bufptr += l;
      break;
    case 'f':
      l=snprintf(bufptr, REMAINING, "f:%g", *(float*)s); s+=4; bufptr += l;
      break;
    case 'd':
      l=snprintf(bufptr, REMAINING, "d:%lg", *(double*)s); s+=8; bufptr += l;
      break;
    case 'z':
    case 'H':
      l=snprintf(bufptr, REMAINING, "%c:", type); bufptr += l;
      while (*s) { 
	if (REMAINING > 0) { *bufptr = *s; bufptr++; }
	s++;
      }
      break;
    case 'B':
      sub_type = *(s++);
      memcpy(&n, s, 4);
      s += 4;
      l=snprintf(bufptr, REMAINING, "%c:%c", type, sub_type); bufptr += l;
      for (j = 0; j < n; ++j) {
	l=snprintf(bufptr, REMAINING, ","); bufptr += l;
	switch (sub_type) {
	case 'c':
	  l=snprintf(bufptr, REMAINING, "%i", *(int8_t*)s); ++s; bufptr += l;
	  break;
	case 'C': 
	  l=snprintf(bufptr, REMAINING, "%u", *(uint8_t*)s); ++s; bufptr += l;
	  break;
	case 's': 
	  l=snprintf(bufptr, REMAINING, "%i", *(int16_t*)s); 
	  s += 2; bufptr += l;
	  break;
	case 'S':
	  l=snprintf(bufptr, REMAINING, "%u", *(uint16_t*)s); 
	  s += 2; bufptr += l;
	  break;
	case 'i': 
	  l=snprintf(bufptr, REMAINING, "%i", *(int32_t*)s); 
	  s += 4; bufptr += l;
	  break;
	case 'I': 
	  l=snprintf(bufptr, REMAINING, "%u", *(uint32_t*)s); 
	  s += 4; bufptr += l;
	  break;
	case 'f':
	  l=snprintf(bufptr, REMAINING, "%g", *(float*)s); 
	  s += 4; bufptr += l;
	  break;
	}
      }
    }
    
    i++;
  }
  SPLICING_CHECK(splicing_strvector_append2(&reads->attributes, buffer, 
					    bufptr-buffer));
  
  /* Single-end or paired-end? */
  
  if (read->core.mpos < 0) {
    reads->noSingles  += 1;
  } else {
    reads->noPairs += 1;
  }
  
  return 0;
}

#undef REMAINING

int splicing_i_order_reads(splicing_reads_t *reads) {
  splicing_vector_char_t taken;
  splicing_vector_int_t idx;
  splicing_vector_t idx2;
  size_t i, pos, noReads=splicing_vector_int_size(&reads->position);
  splicing_i_read_sambam_sort_data_t data = 
    { &reads->chr, &reads->position, &reads->flags, &reads->pairpos };
  int nonpaired_mask = ~0x1;

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
      VECTOR(reads->mypair)[pos] = -1;
      VECTOR(idx2)[pos++] = curr;
    } else if (VECTOR(reads->pairpos)[curr] < VECTOR(reads->position)[curr] 
	       && !VECTOR(taken)[curr])  {
      /* pair is missing */
      VECTOR(reads->mypair)[pos] = -1;
      VECTOR(idx2)[pos++] = curr;
      VECTOR(reads->pairpos)[curr] = 0;
      VECTOR(reads->flags)[curr] &= nonpaired_mask;
      reads->noPairs--;
      reads->noSingles++;
    } else if (VECTOR(reads->pairpos)[curr] < 
	       VECTOR(reads->position)[curr]) {
      /* already done, do nothing */
    } else {
      /* search for pair, forward */
      int needle=VECTOR(reads->pairpos)[curr], ppos, found;
      for (ppos=i+1, found=0; ! found && ppos < noReads; ppos++) {
	int rp=VECTOR(idx)[ppos];
	found = VECTOR(reads->position)[rp] == needle && 
	  VECTOR(reads->pairpos)[rp] == VECTOR(reads->position)[curr] &&
	  0x80 & VECTOR(reads->flags)[rp] && ! VECTOR(taken)[rp];
      }
      if (found) { 
	int rp=VECTOR(idx)[ppos-1];
	VECTOR(reads->mypair)[pos] = pos+1;
	VECTOR(idx2)[pos++] = curr;
	VECTOR(reads->mypair)[pos] = pos-1;
	VECTOR(idx2)[pos++] = rp;
	VECTOR(taken)[rp] = 1;
      } else {
	/* pair is missing */
	VECTOR(reads->mypair)[pos] = -1;
	VECTOR(idx2)[pos++] = curr;
	VECTOR(reads->pairpos)[curr] = 0;
	VECTOR(reads->flags)[curr] &= nonpaired_mask;
	reads->noPairs--;
	reads->noSingles++;
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
  splicing_strvector_permute(&reads->attributes, &idx2);
  
  splicing_vector_destroy(&idx2);
  splicing_vector_int_destroy(&idx);
  SPLICING_FINALLY_CLEAN(2);

  return 0;
}

/* TODO: use qname as well for pairing */

int splicing_read_sambam(const char *filename,
			 splicing_sambam_type_t filetype,
			 splicing_reads_t *reads) {

  samfile_t *infile=0;
  int bytesread;
  bam1_t *read = bam_init1();
  int i;
  char *mode_r="r", *mode_rb="rb", *mode=mode_rb;

  bam_verbose=0;

  SPLICING_FINALLY(splicing_i_bam_destroy1, read);

  splicing_reads_clear(reads);

  switch (filetype) {
    int flen;
  case SPLICING_SAMBAM_AUTO:
    flen=strlen(filename);
    if (flen >= 4 && !strncmp(filename + flen - 4, ".sam", 4)) { 
      mode = mode_r; 
      filetype=SPLICING_SAMBAM_SAM;
    } else {
      filetype=SPLICING_SAMBAM_BAM;
    }
    break;
  case SPLICING_SAMBAM_SAM:
    mode = mode_r;
    break;
  case SPLICING_SAMBAM_BAM:
    mode = mode_rb;
    break;
  }

  infile=samopen(filename, mode, /*aux=*/ 0);
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
    SPLICING_CHECK(splicing_i_add_read(reads, read));
  }

  if (bytesread < -1) {
    SPLICING_WARNING("Truncated SAM/BAM file");
  }

  /* We need to find the pairs of the paired-end reads, if there is any.
     For this we first order the reads according to their RNAME
     (chromosome), plus their position. Then we search for the second
     pair of each first pair. */

  if (reads->noPairs != 0) {
    SPLICING_CHECK(splicing_i_order_reads(reads));
  }

  reads->noPairs /= 2;
  reads->paired = reads->noPairs ? 1 : 0;
    
  bam_destroy1(read);
  SPLICING_FINALLY_CLEAN(1);

  samclose(infile);
  
  return 0;
}

int splicing_i_read_sambam_cb(const bam1_t *read, void *data) {
  splicing_reads_t *reads = (splicing_reads_t*) data;
  SPLICING_CHECK(splicing_i_add_read(reads, read));
}

int splicing_read_sambam_region(const char *filename,
				const char *indexfile,
				splicing_sambam_type_t filetype,
				const char *region,
				splicing_reads_t *reads) {

  samfile_t *infile=0;
  int bytesread;
  bam1_t *read = bam_init1();
  int i;
  char *mode_r="r", *mode_rb="rb", *mode=mode_rb;
  char *myindexfile=(char*) indexfile;
  bam_index_t *idx=0;
  int tid, beg, end, result;

  bam_verbose=0;

  SPLICING_FINALLY(splicing_i_bam_destroy1, read);

  splicing_reads_clear(reads);

  switch (filetype) {
    int flen;
  case SPLICING_SAMBAM_AUTO:
    flen=strlen(filename);
    if (flen >= 4 && !strncmp(filename + flen - 4, ".sam", 4)) { 
      mode = mode_r; 
      filetype=SPLICING_SAMBAM_SAM;
    } else {
      filetype=SPLICING_SAMBAM_BAM;
    }
    break;
  case SPLICING_SAMBAM_SAM:
    mode = mode_r;
    break;
  case SPLICING_SAMBAM_BAM:
    mode = mode_rb;
    break;
  }

  /* Load index, if available */

  if (filetype == SPLICING_SAMBAM_BAM) { 
    FILE * ifp;
    if (!indexfile) { 
      int flen=strlen(filename);
      myindexfile = malloc(flen + 5);
      strcpy(myindexfile, filename);
      strcat(myindexfile, ".bai");
    }
    ifp = fopen(myindexfile, "rb");
    if (!indexfile) { free(myindexfile); }
    if (ifp) { 
      idx = bam_index_load_core(ifp);
      SPLICING_FINALLY(bam_index_destroy, idx);
      fclose(ifp);
    }
    if (!idx) {
      SPLICING_ERROR("Cannot read BAM index file", SPLICING_EFILE);
    }
  }

  infile=samopen(filename, mode, /*aux=*/ 0);
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

  bam_parse_region(infile->header, region, &tid, &beg, &end);
  if (tid < 0) { 
    SPLICING_ERROR("Unknown reference name", SPLICING_EINVAL);
  }
  
  result = bam_fetch(infile->x.bam, idx, tid, beg, end, reads, 
		     splicing_i_read_sambam_cb);
  if (result < 0) {
    SPLICING_ERROR("Truncated BAM file or corrupt BAM index file",
		   SPLICING_EINVAL);
  }

  /* We need to find the pairs of the paired-end reads, if there is any.
     For this we first order the reads according to their RNAME
     (chromosome), plus their position. Then we search for the second
     pair of each first pair. */

  if (reads->noPairs != 0) {
    SPLICING_CHECK(splicing_i_order_reads(reads));
  }

  reads->noPairs /= 2;
  reads->paired = reads->noPairs ? 1 : 0;
  
  if (idx) {
    bam_index_destroy(idx);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  bam_destroy1(read);
  SPLICING_FINALLY_CLEAN(1);

  samclose(infile);  

  return 0;
}
