
#include "splicing_error.h"
#include "splicing.h"

#include <string.h>
#include <ctype.h>
#include <stdio.h>

const char *splicing_types[] = { "gene", "mRNA", "exon", "CDS", 
				 "start_codon", "stop_codon" };

/* TODO: do not ignore size */
int splicing_gff_init(splicing_gff_t *gff, size_t size) {

  if (size < 0) { 
    SPLICING_ERROR("Cannot create GFF, `size' must be non-negative", 
		   SPLICING_EINVAL);
  }

  SPLICING_CHECK(splicing_strvector_init(&gff->seqids, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &gff->seqids);
  SPLICING_CHECK(splicing_strvector_init(&gff->sources, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &gff->sources);
  SPLICING_CHECK(splicing_strvector_init(&gff->ID, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &gff->ID);

  SPLICING_CHECK(splicing_vector_int_init(&gff->genes, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->genes);
  SPLICING_CHECK(splicing_vector_int_init(&gff->transcripts, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->transcripts);
  SPLICING_CHECK(splicing_vector_int_init(&gff->seqid, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->seqid);
  SPLICING_CHECK(splicing_vector_int_init(&gff->source, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->source);
  SPLICING_CHECK(splicing_vector_int_init(&gff->strand, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->strand);
  SPLICING_CHECK(splicing_vector_int_init(&gff->type, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->type);
  SPLICING_CHECK(splicing_vector_int_init(&gff->start, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->start);
  SPLICING_CHECK(splicing_vector_int_init(&gff->end, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->end);
  SPLICING_CHECK(splicing_vector_init(&gff->score, 0));
  SPLICING_FINALLY(splicing_vector_destroy, &gff->score);
  SPLICING_CHECK(splicing_vector_int_init(&gff->phase, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->phase);
  SPLICING_CHECK(splicing_vector_int_init(&gff->parent, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gff->parent);

  gff->n=0;
  gff->nogenes=0;
  gff->notranscripts=0;

  gff->last_gene_id = gff->last_mrna_id = SPLICING_STRVECTOR_ZERO;
  gff->last_gene_no = gff->last_mrna_no = -1;
  gff->last_seqid = gff->last_source = SPLICING_STRVECTOR_ZERO;

  SPLICING_FINALLY_CLEAN(14);

  return 0;
}

void splicing_gff_destroy(splicing_gff_t *gff) {
  splicing_strvector_destroy(&gff->seqids);
  splicing_strvector_destroy(&gff->sources);
  splicing_strvector_destroy(&gff->ID);
  splicing_vector_int_destroy(&gff->genes);
  splicing_vector_int_destroy(&gff->transcripts);
  splicing_vector_int_destroy(&gff->seqid);
  splicing_vector_int_destroy(&gff->source);
  splicing_vector_int_destroy(&gff->strand);
  splicing_vector_int_destroy(&gff->type);
  splicing_vector_int_destroy(&gff->start);
  splicing_vector_int_destroy(&gff->end);
  splicing_vector_destroy(&gff->score);
  splicing_vector_int_destroy(&gff->phase);
  splicing_vector_int_destroy(&gff->parent);
}

void splicing_gff_destroy2(void *gff) {
  splicing_gff_t *o = (splicing_gff_t *) gff;
  splicing_gff_destroy(o);
  free(o);
}

int splicing_gff_reserve(splicing_gff_t *gff, size_t size) {

  SPLICING_CHECK(splicing_vector_int_reserve(&gff->type, size));
  SPLICING_CHECK(splicing_vector_int_reserve(&gff->start, size));
  SPLICING_CHECK(splicing_vector_int_reserve(&gff->end, size));
  SPLICING_CHECK(splicing_vector_reserve(&gff->score, size));
  SPLICING_CHECK(splicing_vector_int_reserve(&gff->phase, size));
  SPLICING_CHECK(splicing_vector_int_reserve(&gff->parent, size));

  return 0;
}

int splicing_gff_append(splicing_gff_t *gff, const char *seqid, 
			const char *source, splicing_type_t type, int start,
			int end, double score, splicing_strand_t strand,
			int phase, const char *ID, const char *parent) {

  if (type == SPLICING_TYPE_GENE) { 
    gff->nogenes++; 
    SPLICING_CHECK(splicing_vector_int_push_back(&gff->genes, gff->n));
    SPLICING_CHECK(splicing_vector_int_push_back(&gff->strand, strand));
  } else if (type == SPLICING_TYPE_MRNA) { 
    gff->notranscripts++; 
    SPLICING_CHECK(splicing_vector_int_push_back(&gff->transcripts, gff->n));
  }

  if (type == SPLICING_TYPE_GENE) {

    /* Seqid */
    if (!strcmp(seqid, gff->last_seqid)) {
      int last=splicing_vector_int_tail(&gff->seqid);
      SPLICING_CHECK(splicing_vector_int_push_back(&gff->seqid, last));
    } else {
      size_t idx;
      int seen=splicing_strvector_search(&gff->seqids, seqid, &idx);
      if (seen) { 
	SPLICING_CHECK(splicing_vector_int_push_back(&gff->seqid, idx));
	gff->last_seqid=splicing_strvector_get(&gff->seqids, idx);
      } else {
	size_t size=splicing_strvector_size(&gff->seqids);
	SPLICING_CHECK(splicing_strvector_append(&gff->seqids, seqid));
	SPLICING_CHECK(splicing_vector_int_push_back(&gff->seqid, size));
	gff->last_source=splicing_strvector_get(&gff->seqids, size);
      }
    }

    /* Source */
    if (!strcmp(source, gff->last_source)) {
      int last=splicing_vector_int_tail(&gff->source);
      SPLICING_CHECK(splicing_vector_int_push_back(&gff->source, last));
    } else {
      size_t idx;
      int seen=splicing_strvector_search(&gff->sources, source, &idx);
      if (seen) { 
	SPLICING_CHECK(splicing_vector_int_push_back(&gff->source, idx));
	gff->last_source=splicing_strvector_get(&gff->sources, idx);
      } else {
	size_t size=splicing_strvector_size(&gff->sources);
	SPLICING_CHECK(splicing_strvector_append(&gff->sources, source));
	SPLICING_CHECK(splicing_vector_int_push_back(&gff->source, size));
	gff->last_source=splicing_strvector_get(&gff->sources, size);
      }
    }

  }

  /* Parent */
  if (!parent || !parent[0]) {
    SPLICING_CHECK(splicing_vector_int_push_back(&gff->parent, -1));
  } else if (!strcmp(parent, gff->last_gene_id)) {
    SPLICING_CHECK(splicing_vector_int_push_back(&gff->parent, 
						 gff->last_gene_no));
  } else if (!strcmp(parent, gff->last_mrna_id)) {
    SPLICING_CHECK(splicing_vector_int_push_back(&gff->parent, 
						 gff->last_mrna_no));
  } else {
    size_t idx;
    int seen=splicing_strvector_search(&gff->ID, parent, &idx);
    if (!seen) { 
      SPLICING_WARNING("Unknown parent ID, invalid GFF file");
      SPLICING_CHECK(splicing_vector_int_push_back(&gff->parent, -1));
    } else {
      SPLICING_CHECK(splicing_vector_int_push_back(&gff->parent, idx));
    }
  }

  SPLICING_CHECK(splicing_vector_int_push_back(&gff->type, type));
  SPLICING_CHECK(splicing_vector_int_push_back(&gff->start, start));
  SPLICING_CHECK(splicing_vector_int_push_back(&gff->end, end));
  SPLICING_CHECK(splicing_vector_push_back(&gff->score, score));
  SPLICING_CHECK(splicing_vector_int_push_back(&gff->phase, phase));
  SPLICING_CHECK(splicing_strvector_append(&gff->ID, ID));
  
  /* Update last gene/mrna */
  if (type == SPLICING_TYPE_GENE) { 
    gff->last_gene_id = splicing_strvector_get(&gff->ID, gff->n);
    gff->last_gene_no = gff->n;
  } else if (type == SPLICING_TYPE_MRNA) {
    gff->last_mrna_id = splicing_strvector_get(&gff->ID, gff->n);
    gff->last_mrna_no = gff->n;
  }

  gff->n += 1;

  return 0;
}

/* Read a string, up to and including the specified delimiter, and put
   it into the given buffer. A zero byte is added to the string.
*/

int splicing_io_get_string(FILE *input, char *buffer, size_t maxlen, 
			   size_t *len, char delim, int newline) {
  int c;

  *len = 0;
  while (1) {
    c=fgetc(input);
    if (c==EOF) {
      SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
    } else if (*len == maxlen) { 
      SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
    } else if (c==delim) {
      *buffer='\0'; 
      buffer++;
      return 0;
    } else if (newline && (c=='\n' || c=='\r')) {
      *buffer='\0';
      buffer++;
      return 0;
    } else if (!newline && (c=='\n' || c=='\r')) {
      SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
    } else { 
      *buffer=(char) c;
      buffer++;
      *len += 1;
    }
  }

  return 1;
}

/* Read an integer, up to and including the specified delimiter. 
   It is an error to have another delimiter character in the input.
 */

int splicing_io_get_integer(FILE *input, int *integer, char delim) {
  char buffer[30];
  char *bufend;
  size_t len;
  int eof = splicing_io_get_string(input, buffer, 
				   sizeof(buffer)/sizeof(char), &len, 
				   delim, /*newline=*/ 0);

  if (eof) { 
    SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
  }

  *integer = (int) strtol(buffer, &bufend, /*base=*/ 10);
  if (*bufend != '\0') { 
    SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
  }

  return 0;
}

int splicing_io_get_integer_na(FILE *input, int *integer, char delim,
			       char nachar) {
  char buffer[30];
  char *bufend;
  size_t len;
  int na=SPLICING_NA_INTEGER;
  int eof = splicing_io_get_string(input, buffer, 
				   sizeof(buffer)/sizeof(char), &len, 
				   delim, /*newline=*/ 0);

  if (eof) { 
    SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
  }

  if (len > 0 && buffer[0]==nachar) {
    *integer=na;
    return 0;
  }

  *integer = (int) strtol(buffer, &bufend, /*base=*/ 10);
  if (*bufend != '\0') { 
    SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
  }
  
  return 0;
}

/* Read a real value, up to and including the specified delimiter. 
   It is an error to have another delimiter character in the input.
 */

int splicing_io_get_real_na(FILE *input, double *real, char delim, 
			    char nachar) {
  char buffer[30];
  char *bufend;
  size_t len;
  double na=SPLICING_NA_REAL;
  int eof = splicing_io_get_string(input, buffer, 
				   sizeof(buffer)/sizeof(char), &len, 
				   delim, /*newline=*/ 0);

  if (eof) { 
    SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
  }

  if (len > 0 && buffer[0]==nachar) {
    *real=na;
    return 0;
  }

  *real = strtod(buffer, &bufend);
  if (*bufend != '\0') { 
    SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
  }

  return 0;
}

/* Skip the newline character(s) on the input, 
   plus any lines starting with a '#' character */

int splicing_io_skip_newline_and_comments(FILE *input) {
  int c;

  do {
    c=fgetc(input);
    if (c=='#') { 
      do {
	c=fgetc(input);
      } while (c != '\n' && c != '\r');
    }
  } while (c != EOF && (c=='\n' || c=='\r'));

  if (c != EOF) { 
    ungetc(c, input);
    return 0;
  } else {
    return EOF;
  }
}

int splicing_io_parse_attributes(char *attr, char **ID, char**parent) {
  *ID=SPLICING_STRVECTOR_ZERO; *parent=SPLICING_STRVECTOR_ZERO;
  char *kw, *vl;
  while (*attr != '\0') {
    /* Skip white space */
    while (*attr != '\0' && isspace(*attr)) {
      attr++;
    }
    /* Keyword */
    kw=attr;
    while (*attr != '\0' && *attr != '=') {
      attr++;
    }
    if (*attr == '\0') { 
      SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
    }
    *attr='\0'; attr++;
    /* Value */
    vl=attr;
    while (*attr != '\0' && *attr != ';') { 
      attr++;
    }
    if (*attr == ';') { *attr='\0'; attr++; }
    if (!strcmp("ID", kw)) {
      *ID=vl;
    } else if (!strcmp("Parent", kw)) {
      *parent=vl;
    }
  }
  return 0;
}

/* TODO: sort it */

int splicing_gff_read(FILE *input, splicing_gff_t *gff) {
  int eof=!EOF;
  char seqid[200];
  char source[200];
  char type[200];
  int start, end, phase;
  double score;
  char strand[10];
  char attributes[5000];
  size_t len;
  splicing_type_t realtype;
  splicing_strand_t realstrand=SPLICING_STRAND_UNKNOWN;
  char *ID, *parent;

  do { 

    eof = eof || splicing_io_skip_newline_and_comments(input);
    eof = eof || splicing_io_get_string(input, seqid, 
					sizeof(seqid)/sizeof(char), &len, 
					/*delim=*/ '\t', /*newline=*/ 0);
    eof = eof || splicing_io_get_string(input, source, 
					sizeof(source)/sizeof(char), &len,
					/*delim=*/ '\t', /*newline=*/ 0);
    eof = eof || splicing_io_get_string(input, type, 
					sizeof(type)/sizeof(char), &len, 
					/*delim=*/ '\t', /*newline=*/ 0);
    eof = eof || splicing_io_get_integer(input, &start, /*delim=*/ '\t');
    eof = eof || splicing_io_get_integer(input, &end, /*delim=*/ '\t');
    eof = eof || splicing_io_get_real_na(input, &score, /*delim=*/ '\t', 
					 /*nachar=*/ '.');
    eof = eof || splicing_io_get_string(input, strand, 
					sizeof(strand)/sizeof(char), &len, 
					/*delim=*/ '\t', /*newline=*/ 0);
    eof = eof || splicing_io_get_integer_na(input, &phase, /*delim=*/ '\t', 
					    /*nachar=*/ '.');
    eof = eof || splicing_io_get_string(input, attributes, 
					sizeof(attributes)/sizeof(char), &len,
					/*delim=*/ '\n', /*newline=*/ 1);

    if (eof) { 
      SPLICING_ERROR("Corrupt GFF file", SPLICING_PARSEERROR);
    }

    /* TODO: do not hardcode these names
       TODO: order them according to their frequency */
    if (!strcmp(type, "gene")) { 
      realtype = SPLICING_TYPE_GENE;
    } else if (!strcmp(type, "mRNA")) {
      realtype = SPLICING_TYPE_MRNA;
    } else if (!strcmp(type, "exon")) {
      realtype = SPLICING_TYPE_EXON;
    } else if (!strcmp(type, "CDS")) {
      realtype = SPLICING_TYPE_CDS;
    } else if (!strcmp(type, "start_codon")) {
      realtype = SPLICING_TYPE_START_CODON;
    } else if (!strcmp(type, "stop_codon")) {
      realtype = SPLICING_TYPE_STOP_CODON;
    } else {
      SPLICING_ERROR("Invalid GFF file", SPLICING_PARSEERROR);
    }

    if (!strcmp(strand, "+")) {
      realstrand=SPLICING_STRAND_PLUS;
    } else if (!strcmp(strand, "-")) {
      realstrand=SPLICING_STRAND_MINUS;
    } else if (!strcmp(strand, ".")) {
      realstrand=SPLICING_STRAND_UNKNOWN;
    }

    /* Parsing the attributes field */
    SPLICING_CHECK(splicing_io_parse_attributes(attributes, &ID, &parent));

    SPLICING_CHECK(splicing_gff_append(gff, seqid, source, realtype, start,
				       end, score, realstrand, phase, ID,
				       parent));
    
    eof = splicing_io_skip_newline_and_comments(input);

  } while (!eof);
  
  return 0;
}

int splicing_gff_nogenes(const splicing_gff_t *gff, size_t *nogenes) {
  *nogenes = splicing_vector_int_size(&gff->genes);
  return 0;
}

int splicing_i_gff_noiso_one(const splicing_gff_t *gff, size_t gene,
			     size_t *noiso, splicing_vector_int_t *isolen) {

  size_t nogenes, idx1, idx2;
  SPLICING_CHECK(splicing_gff_nogenes(gff, &nogenes));
  
  if (gene < 0 || gene >= nogenes) {
    SPLICING_ERROR("Invalid gene id", SPLICING_EINVAL);
  }

  idx1=VECTOR(gff->genes)[gene];
  idx2= gene+1 == nogenes ? gff->n : VECTOR(gff->genes)[gene+1];
  
  *noiso = 0;
  for ( ; idx1 < idx2; idx1++) {
    if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_MRNA) { *noiso += 1; }
  }

  if (isolen) {
    size_t il=0, pos=0;
    SPLICING_CHECK(splicing_vector_int_resize(isolen, *noiso));
    idx1=VECTOR(gff->genes)[gene];
    idx2= gene+1 == nogenes ? gff->n : VECTOR(gff->genes)[gene+1];
    
    for (; idx1 < idx2 && VECTOR(gff->type)[idx1] != SPLICING_TYPE_MRNA; 
	 idx1++) ;
    idx1++;
    for (; idx1 < idx2; idx1++) {
      if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_MRNA) { 
	VECTOR(*isolen)[pos++]=il;
	il = 0;
      } else if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_EXON) {
	il += VECTOR(gff->end)[idx1] - VECTOR(gff->start)[idx1] + 1;
      }
    }
    VECTOR(*isolen)[pos++]=il;
  }
  
  return 0;
}

int splicing_gff_noexons_one(const splicing_gff_t *gff, size_t gene,
			     splicing_vector_int_t *noexons) {

  size_t nogenes, idx1, idx2, noiso, pos, il;
  SPLICING_CHECK(splicing_gff_nogenes(gff, &nogenes));
  
  if (gene < 0 || gene >= nogenes) {
    SPLICING_ERROR("Invalid gene id", SPLICING_EINVAL);
  }

  idx1=VECTOR(gff->genes)[gene];
  idx2= gene+1 == nogenes ? gff->n : VECTOR(gff->genes)[gene+1];
  
  for (noiso=0; idx1 < idx2; idx1++) {
    if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_MRNA) { noiso += 1; }    
  }

  SPLICING_CHECK(splicing_vector_int_resize(noexons, noiso));

  idx1=VECTOR(gff->genes)[gene];
  idx2= gene+1 == nogenes ? gff->n : VECTOR(gff->genes)[gene+1];
  for (; idx1 < idx2 && VECTOR(gff->type)[idx1] != SPLICING_TYPE_MRNA; 
       idx1++) ;
  idx1++;
  for (pos=0, il=0; idx1 < idx2; idx1++) {
    if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_MRNA) {
      VECTOR(*noexons)[pos++]=il;
      il=0;
    } else if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_EXON) { il++; }
  }
  VECTOR(*noexons)[pos++]=il;

  return 0;
}

/* TODO: error checks */

int splicing_gff_noiso(const splicing_gff_t *gff, 
		       splicing_vector_int_t *noiso) {

  size_t nogenes, idx1, idx2, pos=0;

  SPLICING_CHECK(splicing_gff_nogenes(gff, &nogenes));
  idx1=VECTOR(gff->genes)[0];
  idx2=gff->n;
  
  SPLICING_CHECK(splicing_vector_int_resize(noiso, nogenes));
  splicing_vector_int_null(noiso);
  if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_GENE) { idx1++; }
  for (; idx1 < gff->n; idx1++) { 
    if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_MRNA) { 
      VECTOR(*noiso)[pos] += 1;
    } else if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_GENE) {
      pos++;
    }
  }

  return 0;
}

int splicing_gff_noiso_one(const splicing_gff_t *gff, size_t gene, 
			   size_t *noiso) {
  return splicing_i_gff_noiso_one(gff, gene, noiso, 0);
}

int splicing_gff_isolength_one(const splicing_gff_t *gff, size_t gene, 
			       splicing_vector_int_t *isolength) {
  size_t noiso;
  return splicing_i_gff_noiso_one(gff, gene, &noiso, isolength);
}

int splicing_gff_isolength(const splicing_gff_t *gff,
			   splicing_vector_int_t *isolength,
			   splicing_vector_int_t *isolength_idx) {

  size_t idx1;
  size_t nogenes=splicing_vector_int_size(&gff->genes);
  size_t notrans=splicing_vector_int_size(&gff->transcripts);
  int pos=-1, ipos=0;
  
  SPLICING_CHECK(splicing_vector_int_resize(isolength, notrans));
  SPLICING_CHECK(splicing_vector_int_resize(isolength_idx, nogenes));
  
  for (idx1=VECTOR(gff->genes)[0]; idx1 < gff->n; idx1++) {
    if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_GENE) {
      VECTOR(*isolength_idx)[ipos++]=pos+1;
    } else if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_MRNA) {
      VECTOR(*isolength)[++pos] = 0;
    } else if (VECTOR(gff->type)[idx1] == SPLICING_TYPE_EXON) {
      VECTOR(*isolength)[pos] +=
	(VECTOR(gff->end)[idx1] - VECTOR(gff->start)[idx1] + 1);
    }
  }

  return 0;
}

size_t splicing_gff_size(const splicing_gff_t *gff) {
  return gff->n;
}

int splicing_i_gff_exon_start_end_sort_cmp(void *data, const void *a, 
					   const void *b) {

  int *start=(int*) data;
  int aa=*(int*)a, bb=*(int*)b;
  int sa=start[aa];
  int sb=start[bb];
  
  if (sa < sb) {
    return -1;
  } else if (sa > sb) { 
    return 1;
  } else {
    return 0;
  }
}

int splicing_i_gff_exon_start_end_sort(const splicing_vector_int_t *start,
				       const splicing_vector_int_t *end,
				       const splicing_vector_int_t *idx, 
				       int iso, splicing_vector_int_t *tmp, 
				       splicing_vector_int_t *tmp2) {
  
  int i, j, from=VECTOR(*idx)[iso], to=VECTOR(*idx)[iso+1], len=to-from;  

  SPLICING_CHECK(splicing_vector_int_resize(tmp, len));
  SPLICING_CHECK(splicing_vector_int_resize(tmp2, len));
  for (i=0; i<len; i++) { VECTOR(*tmp)[i]=i; }
  splicing_qsort_r(VECTOR(*tmp), len, sizeof(int), 
		   (void*) (VECTOR(*start)+from),
		   splicing_i_gff_exon_start_end_sort_cmp);
  
  /* Store the order */
  for (i=0, j=from; i<len; i++, j++) { VECTOR(*tmp2)[i]=VECTOR(*start)[j]; }
  for (i=0, j=from; i<len; i++, j++) { 
    VECTOR(*start)[j] = VECTOR(*tmp2)[ VECTOR(*tmp)[i] ];
  }
  for (i=0, j=from; i<len; i++, j++) { VECTOR(*tmp2)[i]=VECTOR(*end)[j]; }
  for (i=0, j=from; i<len; i++, j++) { 
    VECTOR(*end)[j] = VECTOR(*tmp2)[ VECTOR(*tmp)[i] ];
  }

  return 0;
}

/* Return the start coordinates of all exons, in each isoform, 
   for a given gene. This function also makes sure that the 
   exons are ordered according to their start coordinate, within each
   isoform. */

int splicing_gff_exon_start_end(const splicing_gff_t *gff, 
				splicing_vector_int_t *start,
				splicing_vector_int_t *end,
				splicing_vector_int_t *idx,
				int gene) {
  
  size_t noiso;
  int i=0, p=0, n=splicing_gff_size(gff);
  int pos;
  size_t nogenes;
  splicing_vector_int_t tmp, tmp2;

  SPLICING_CHECK(splicing_vector_int_init(&tmp, 10));
  SPLICING_FINALLY(splicing_vector_int_destroy, &tmp);
  SPLICING_CHECK(splicing_vector_int_init(&tmp2, 10));
  SPLICING_FINALLY(splicing_vector_int_destroy, &tmp2);

  SPLICING_CHECK(splicing_gff_nogenes(gff, &nogenes));
  if (gene < 0 || gene >= nogenes) { 
    SPLICING_ERROR("Invalid gene id", SPLICING_EINVAL);
  }

  pos=VECTOR(gff->genes)[gene]+1;
  
  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));
  splicing_vector_int_clear(start);
  splicing_vector_int_clear(end);
  SPLICING_CHECK(splicing_vector_int_resize(idx, noiso+1));
  while (pos < n) {
    if (VECTOR(gff->type)[pos] == SPLICING_TYPE_EXON) { 
      int s=VECTOR(gff->start)[pos];
      int e=VECTOR(gff->end)[pos];
      SPLICING_CHECK(splicing_vector_int_push_back(start, s)); p++;
      SPLICING_CHECK(splicing_vector_int_push_back(end, e));
    } else if (VECTOR(gff->type)[pos] == SPLICING_TYPE_MRNA) {
      VECTOR(*idx)[i] = p;
      if (i!=0) { 
	SPLICING_CHECK(splicing_i_gff_exon_start_end_sort(start, end, idx, 
							  i-1, &tmp, &tmp2));
      }
      i++;
    } else if (VECTOR(gff->type)[pos] == SPLICING_TYPE_GENE) {
      break;
    }
    pos++;
  }
  VECTOR(*idx)[i] = p;
  SPLICING_CHECK(splicing_i_gff_exon_start_end_sort(start, end, idx, i-1, 
						    &tmp, &tmp2));

  splicing_vector_int_destroy(&tmp2);
  splicing_vector_int_destroy(&tmp);
  SPLICING_FINALLY_CLEAN(2);

  return 0;
}

/* TODO: handle errors */

int splicing_gff_write(FILE *output, const splicing_gff_t *gff) {

  size_t i, n=gff->n;
  size_t gene=-1;
  const char *types[] = { "gene", "mRNA", "exon", "CDS", "start_codon",
			  "stop_codon" };
  const char *strands[] = { "+", "-", "." };
  
  for (i=0; i<n; i++) {

    if (VECTOR(gff->type)[i] == SPLICING_TYPE_GENE) { gene++; }

    fprintf(output, "%s\t%s\t%s\t%li\t%li\t", 
	    splicing_strvector_get(&gff->seqids, VECTOR(gff->seqid)[gene]),
	    splicing_strvector_get(&gff->sources, VECTOR(gff->source)[gene]),
	    types[VECTOR(gff->type)[i]], (long) VECTOR(gff->start)[i],
	    (long) VECTOR(gff->end)[i]);
    if (VECTOR(gff->score)[i] == SPLICING_NA_REAL) { 
      fprintf(output, "%s\t", ".");
    } else {
      fprintf(output, "%g\t", VECTOR(gff->score)[i]);
    }
    fprintf(output, "%s\t", strands[gene]);
    if (VECTOR(gff->phase)[i] == SPLICING_NA_INTEGER) {
      fprintf(output, "%s\t", ".");
    } else {
      fprintf(output, "%i\t", VECTOR(gff->phase)[i]);
    }
    fprintf(output, "ID=%s", splicing_strvector_get(&gff->ID, i));
    if (VECTOR(gff->parent)[i] >= 0) {
      fprintf(output, ";Parent=%s", 
	      splicing_strvector_get(&gff->ID, VECTOR(gff->parent)[i]));
    }
    fputs("\n", output);
  }
  
  return 0;
}

int splicing_gff_gene_start_end_one(const splicing_gff_t *gff, size_t gene,
				    size_t *start, size_t *end) {

  size_t nogenes=splicing_vector_int_size(&gff->genes);
  size_t idx;
  
  if (gene < 0 || gene >= nogenes) { 
    SPLICING_ERROR("Invalid gene id", SPLICING_EINVAL); 
  }
  
  idx=VECTOR(gff->genes)[gene];
  *start=VECTOR(gff->start)[idx];
  *end=VECTOR(gff->end)[idx];

  return 0;
}

int splicing_gff_gene_start_end(const splicing_gff_t *gff, 
				splicing_vector_int_t *start,
				splicing_vector_int_t *end) {

  size_t i, nogenes=splicing_vector_int_size(&gff->genes);
  
  SPLICING_CHECK(splicing_vector_int_resize(start, nogenes));
  SPLICING_CHECK(splicing_vector_int_resize(end, nogenes));
  
  for (i=0; i<nogenes; i++) {
    size_t idx=VECTOR(gff->genes)[i];
    VECTOR(*start)[i] = VECTOR(gff->start)[idx];
    VECTOR(*end)[i] = VECTOR(gff->end)[idx];
  }
  
  return 0;
}

int splicing_iso_to_genomic(const splicing_gff_t *gff, size_t gene, 
			    const splicing_vector_int_t *isoform,
			    const splicing_vector_int_t *exstart,
			    const splicing_vector_int_t *exend,
			    const splicing_vector_int_t *exidx,
			    splicing_vector_int_t *position) {

  size_t i, noiso, n=splicing_vector_int_size(position);
  splicing_vector_int_t exlim, shift;
  splicing_vector_int_t vexstart, vexend, vexidx, 
    *myexstart=(splicing_vector_int_t *) exstart, 
    *myexend=(splicing_vector_int_t *) exend, 
    *myexidx=(splicing_vector_int_t *) exidx;
  size_t pos, pos2;

  if (!exstart || !exend || !exidx) {
    myexstart=&vexstart;
    myexend=&vexend;
    myexidx=&vexidx;
    SPLICING_CHECK(splicing_vector_int_init(myexstart, 0));
    SPLICING_FINALLY(splicing_vector_int_destroy, myexstart);
    SPLICING_CHECK(splicing_vector_int_init(myexend, 0));
    SPLICING_FINALLY(splicing_vector_int_destroy, myexend);
    SPLICING_CHECK(splicing_vector_int_init(myexidx, 0));
    SPLICING_FINALLY(splicing_vector_int_destroy, myexidx);
    SPLICING_CHECK(splicing_gff_exon_start_end(gff, myexstart, myexend, 
					       myexidx, gene));
  }

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));

  SPLICING_CHECK(splicing_vector_int_init(&exlim, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exlim);
  SPLICING_CHECK(splicing_vector_int_init(&shift, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &shift);

  for (i=0; i<noiso; i++) {
    size_t cs=0, ce=0, ex=0;
    int pos=VECTOR(*myexidx)[i], pos2=VECTOR(*myexidx)[i+1];
    while (pos < pos2) {
      cs += VECTOR(*myexstart)[pos];
      SPLICING_CHECK(splicing_vector_int_push_back(&shift, cs-ce-ex-1));
      ex++; ce += VECTOR(*myexend)[pos]; pos++;
    }
  }

  for (i=0; i<noiso; i++) { 
    size_t cs=0;
    int pos=VECTOR(*myexidx)[i], pos2=VECTOR(*myexidx)[i+1];
    while (pos < pos2) {
      size_t l=VECTOR(*myexend)[pos]-VECTOR(*myexstart)[pos]+1;
      cs += l;
      SPLICING_CHECK(splicing_vector_int_push_back(&exlim, cs+1));
      pos++;
    }
  }  

  for (i=0; i<n; i++) {
    int iso=VECTOR(*isoform)[i];
    size_t pos=VECTOR(*position)[i];
    int ex;
    for (ex=VECTOR(*myexidx)[iso]; VECTOR(exlim)[ex] <= pos; ex++) ;
    VECTOR(*position)[i] = pos + VECTOR(shift)[ex];
  }

  splicing_vector_int_destroy(&shift);
  splicing_vector_int_destroy(&exlim);
  SPLICING_FINALLY_CLEAN(2);

  if (!exstart || !exend || !exidx) {
    splicing_vector_int_destroy(myexidx);
    splicing_vector_int_destroy(myexend);
    splicing_vector_int_destroy(myexstart);
    SPLICING_FINALLY_CLEAN(3);
  }
  
  return 0;
}

int splicing_genomic_to_iso(const splicing_gff_t *gff, size_t gene,
			    const splicing_vector_int_t *position, 
			    splicing_matrix_int_t *isopos) {

  size_t r, i, noiso, noreads=splicing_vector_int_size(position);
  splicing_vector_int_t exstart, exend, exidx, shift;
  
  splicing_gff_noiso_one(gff, gene, &noiso);
  
  SPLICING_CHECK(splicing_vector_int_init(&exstart, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exstart);
  SPLICING_CHECK(splicing_vector_int_init(&exend, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exend);
  SPLICING_CHECK(splicing_vector_int_init(&exidx, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &exidx);
  SPLICING_CHECK(splicing_gff_exon_start_end(gff, &exstart, &exend,
					     &exidx, gene));

  SPLICING_CHECK(splicing_vector_int_init(&shift, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &shift);
  
  for (i=0; i<noiso; i++) {
    size_t cs=0, ce=0, ex=0;
    int pos=VECTOR(exidx)[i], pos2=VECTOR(exidx)[i+1];
    while (pos < pos2) {
      cs += VECTOR(exstart)[pos];
      SPLICING_CHECK(splicing_vector_int_push_back(&shift, cs-ce-ex-1));
      ex++; ce += VECTOR(exend)[pos]; pos++;
    }
  }

  SPLICING_CHECK(splicing_matrix_int_resize(isopos, noiso, noreads));
  
  for (r=0; r<noreads; r++) {
    for (i=0; i<noiso; i++) {
      size_t pos=VECTOR(*position)[r];
      size_t startpos=VECTOR(exidx)[i];
      size_t endpos=VECTOR(exidx)[i+1];
      int ex;
      for (ex=startpos; ex < endpos && VECTOR(exend)[ex] < pos; ex++) ;
      if (VECTOR(exstart)[ex] <= pos && pos <= VECTOR(exend)[ex]) {
	MATRIX(*isopos, i, r) = VECTOR(*position)[r] - VECTOR(shift)[ex];
      } else { 
	MATRIX(*isopos, i, r) = -1;
      }
    }
  }

  splicing_vector_int_destroy(&shift);
  splicing_vector_int_destroy(&exidx);
  splicing_vector_int_destroy(&exend);
  splicing_vector_int_destroy(&exstart);
  SPLICING_FINALLY_CLEAN(4);

  return 0;
}
	
int splicing_gff_fprint_gene(const splicing_gff_t *gff, 
			     FILE *outfile, int gene) {

  size_t nogenes, noiso;
  int i, j;
  splicing_vector_int_t start, end, idx;

  SPLICING_CHECK(splicing_gff_nogenes(gff, &nogenes));
  
  if (gene < 0 || gene >= nogenes) { 
    SPLICING_ERROR("Invalid gene ID", SPLICING_EINVAL);
  }

  SPLICING_CHECK(splicing_vector_int_init(&start, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &start);
  SPLICING_CHECK(splicing_vector_int_init(&end, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &end);
  SPLICING_CHECK(splicing_vector_int_init(&idx, 0));  
  SPLICING_FINALLY(splicing_vector_int_destroy, &idx);

  SPLICING_CHECK(splicing_gff_exon_start_end(gff, &start, &end, &idx, gene));
  noiso = splicing_vector_int_size(&idx)-1;
  
  fprintf(outfile, "===\nGene with %i isoforms:\n", (int) noiso);
  for (i=0; i<noiso; i++) {
    fprintf(outfile, "  Isoform %i:\n", i);
    for (j=VECTOR(idx)[i]; j<VECTOR(idx)[i+1]; j++) {
      fprintf(outfile, "    %i-%i\n", VECTOR(start)[j], VECTOR(end)[j]);
    }
  }
  
  splicing_vector_int_destroy(&idx);
  splicing_vector_int_destroy(&end);
  splicing_vector_int_destroy(&start);
  SPLICING_FINALLY_CLEAN(3);
  
  return 0;    
}

int splicing_gff_print_gene(const splicing_gff_t *gff, 
			    int gene) {
  return splicing_gff_fprint_gene(gff, stdout, gene);
}
		    
int splicing_gff_fprint(const splicing_gff_t *gff, 
			FILE *outfile) {

  size_t i, n;
  SPLICING_CHECK(splicing_gff_nogenes(gff, &n));
  for (i=0; i<n; i++) {
    SPLICING_CHECK(splicing_gff_fprint_gene(gff, outfile, i));
  }
  
  return 0;
}

int splicing_gff_print(const splicing_gff_t *gff) {
  return splicing_gff_fprint(gff, stdout);
}
