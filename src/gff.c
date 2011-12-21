
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

#define STR(x) (splicing_strvector_get(&gff->ID, (x)))

int splicing_i_gff_reindex_cmp(void *data, const void *a, const void *b) {
  splicing_gff_t *gff=(splicing_gff_t *) data;
  int aa=*(int*)a, bb=*(int*)b;
  
  int parent_a=VECTOR(gff->parent)[aa];
  int parent_b=VECTOR(gff->parent)[bb];
  int gparent_a= parent_a == -1 ? -1 : VECTOR(gff->parent)[parent_a];
  int gparent_b= parent_b == -1 ? -1 : VECTOR(gff->parent)[parent_b];
  
  const char *a_gene_id, *b_gene_id, *a_mrna_id, *b_mrna_id;
  int c1, c2;

  /* If gene ids differ */
  a_gene_id = gparent_a != -1 ? STR(gparent_a) : 
    (parent_a != -1 ? STR(parent_a) : STR(aa));
  b_gene_id = gparent_b != -1 ? STR(gparent_b) : 
    (parent_b != -1 ? STR(parent_b) : STR(bb));
  c1=strcmp(a_gene_id, b_gene_id); if (c1 != 0) { return c1; }

  /* Or if mRNA ids differ */
  a_mrna_id = gparent_a != -1 ? STR(parent_a) : STR(aa);
  b_mrna_id = gparent_b != -1 ? STR(parent_b) : STR(bb);
  c2=strcmp(a_mrna_id, b_mrna_id); if (c2 != 0) { return c2; }
  
  /* Otherwise gene first, then mRNA, then the rest according to
     start position */
  if (parent_a == -1 && parent_b != -1) { 
    return -1; 
  } else if (parent_a != -1 && parent_b == -1) { 
    return 1;
  } else if (gparent_a == -1 && gparent_b != -1) { 
    return -1;
  } else if (gparent_a != -1 && gparent_b == -1) { 
    return 1;
  } else if (gparent_a != -1 && gparent_b != -1) { 
    int sa=VECTOR(gff->start)[aa];
    int sb=VECTOR(gff->start)[bb];
    if (sa < sb) { return -1; } else if (sa > sb) { return 1; }
    return 0;
  } else {
    SPLICING_ERROR("Invalid GFF file, cannot order records", 
		   SPLICING_EINVAL);
  }
  return 0;
}

#undef STR

int splicing_gff_reindex(splicing_gff_t *gff) {
  splicing_vector_int_t index;
  splicing_vector_int_t index2;
  splicing_vector_int_t gindex;
  int i, j, k, n=gff->n;
  
  SPLICING_CHECK(splicing_vector_int_init(&index, n));
  SPLICING_FINALLY(splicing_vector_int_destroy, &index);

  for (i=0; i<n; i++) { VECTOR(index)[i] = i; }

  splicing_qsort_r(VECTOR(index), n, sizeof(int), (void*) gff, 
		   splicing_i_gff_reindex_cmp);

  SPLICING_CHECK(splicing_vector_int_init(&gindex, gff->nogenes));
  SPLICING_FINALLY(splicing_vector_int_destroy, &gindex);
  SPLICING_CHECK(splicing_vector_int_init(&index2, n));
  SPLICING_FINALLY(splicing_vector_int_destroy, &index2);

  for (i=0; i<gff->nogenes; i++) {
    VECTOR(index2)[ VECTOR(gff->genes)[i] ] = i;
  }

  for (i=0, j=0; i<n; i++) {
    if (VECTOR(gff->type)[ VECTOR(index)[i] ] == SPLICING_TYPE_GENE) {
      VECTOR(gindex)[j++] = VECTOR(index2)[ VECTOR(index)[i] ];
    }
  }

  splicing_vector_int_destroy(&index2);
  SPLICING_FINALLY_CLEAN(1);

  splicing_vector_int_intiindex(&gff->seqid, &gindex);
  splicing_vector_int_intiindex(&gff->source, &gindex);
  splicing_vector_int_intiindex(&gff->strand, &gindex);  
    
  splicing_vector_int_destroy(&gindex);
  SPLICING_FINALLY_CLEAN(1);

  splicing_vector_int_intiindex(&gff->type, &index);
  splicing_vector_int_intiindex(&gff->start, &index);
  splicing_vector_int_intiindex(&gff->end, &index);
  splicing_vector_intiindex(&gff->score, &index);
  splicing_vector_int_intiindex(&gff->phase, &index);
  splicing_strvector_ipermute(&gff->ID, &index);
  splicing_vector_int_intiindex(&gff->parent, &index);

  SPLICING_CHECK(splicing_vector_int_init(&index2, n));
  SPLICING_FINALLY(splicing_vector_int_destroy, &index2);

  for (i=0; i<n; i++) {
    VECTOR(index2)[ VECTOR(index)[i] ] = i;
  }

  for (i=0; i<n; i++) {
    int p=VECTOR(gff->parent)[i];
    if (p != -1) {
      VECTOR(gff->parent)[i] = VECTOR(index2)[p];
    }
  }

  splicing_vector_int_destroy(&index2);
  SPLICING_FINALLY_CLEAN(1);

  for (i=j=k=0; i<n; i++) { 
    if (VECTOR(gff->type)[i] == SPLICING_TYPE_GENE) { 
      VECTOR(gff->genes)[j++] = i;
    } else if (VECTOR(gff->type)[i] == SPLICING_TYPE_MRNA) { 
      VECTOR(gff->transcripts)[k++] = i;
    }
  }  

  splicing_vector_int_destroy(&index);
  SPLICING_FINALLY_CLEAN(1);

  return 0;
}

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

  SPLICING_CHECK(splicing_gff_reindex(gff));

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

/* Return the start coordinates of all exons, in each isoform, 
   for a given gene. We assume that the gff is correctly sorted. */

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
      i++;
    } else if (VECTOR(gff->type)[pos] == SPLICING_TYPE_GENE) {
      break;
    }
    pos++;
  }
  VECTOR(*idx)[i] = p;

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
    fprintf(output, "%s\t", strands[VECTOR(gff->strand)[gene]]);
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

int splicing_gff_converter_init(const splicing_gff_t *gff, size_t gene,
				splicing_gff_converter_t *converter) {

  int i; 

  SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &converter->noiso));

  SPLICING_VECTOR_INT_INIT_FINALLY(&converter->exstart, 0);
  SPLICING_VECTOR_INT_INIT_FINALLY(&converter->exend, 0);
  SPLICING_VECTOR_INT_INIT_FINALLY(&converter->exidx, 0);
  SPLICING_VECTOR_INT_INIT_FINALLY(&converter->shift, 0);
  SPLICING_VECTOR_INT_INIT_FINALLY(&converter->exlim, 0);
  
  SPLICING_CHECK(splicing_gff_exon_start_end(gff, &converter->exstart, 
					     &converter->exend, 
					     &converter->exidx, gene));

  /* Calculate the shift */
  for (i=0; i < converter->noiso; i++) {
    size_t cs=0, ce=0, ex=0;
    int pos=VECTOR(converter->exidx)[i], pos2=VECTOR(converter->exidx)[i+1];
    while (pos < pos2) {
      cs += VECTOR(converter->exstart)[pos];
      SPLICING_CHECK(splicing_vector_int_push_back(&converter->shift, 
						   cs-ce-ex-1));
      ex++; ce += VECTOR(converter->exend)[pos]; pos++;
    }
  }
  
  /* Calculate the exlim */
  for (i=0; i < converter->noiso; i++) { 
    size_t cs=0;
    int pos=VECTOR(converter->exidx)[i], pos2=VECTOR(converter->exidx)[i+1];
    while (pos < pos2) {
      size_t l=
	VECTOR(converter->exend)[pos] - VECTOR(converter->exstart)[pos]+1;
      cs += l;
      SPLICING_CHECK(splicing_vector_int_push_back(&converter->exlim, cs+1));
      pos++;
    }
  }

  SPLICING_FINALLY_CLEAN(5);

  return 0;
}

void splicing_gff_converter_destroy(splicing_gff_converter_t *converter) {
  splicing_vector_int_destroy(&converter->exstart);
  splicing_vector_int_destroy(&converter->exend);
  splicing_vector_int_destroy(&converter->exidx);
  splicing_vector_int_destroy(&converter->exlim);
  splicing_vector_int_destroy(&converter->shift);
}

/* Convert isoform coordinates to genomic coordinates. 
   Convert a vector of positions (position), this will be overwritten
   by the result. For each position, its isoform can be specified
   (isoform). 

   If input position is negative, then it is ignored.
*/

int splicing_iso_to_genomic(const splicing_gff_t *gff, size_t gene, 
			    const splicing_vector_int_t *isoform,
			    const splicing_gff_converter_t *converter,
			    splicing_vector_int_t *position) {

  size_t i, n=splicing_vector_int_size(position);
  splicing_gff_converter_t vconverter, 
    *myconverter = (splicing_gff_converter_t*) converter;

  if (!converter) { 
    myconverter=&vconverter;
    SPLICING_CHECK(splicing_gff_converter_init(gff, gene, myconverter));
    SPLICING_FINALLY(splicing_gff_converter_destroy, myconverter);
  }

  /* Do the shifting */
  for (i=0; i<n; i++) {
    int iso=VECTOR(*isoform)[i];
    size_t pos=VECTOR(*position)[i];
    int ex;
    if (pos==-1) { continue; }
    for (ex=VECTOR(myconverter->exidx)[iso]; 
	 ex < VECTOR(myconverter->exidx)[iso+1] && 
	   VECTOR(myconverter->exlim)[ex] <= pos; 
	 ex++) ;
    if (ex < VECTOR(myconverter->exidx)[iso+1]) { 
      VECTOR(*position)[i] = pos + VECTOR(myconverter->shift)[ex];
    } else {
      VECTOR(*position)[i] = -1;
    }
  }

  if (!converter) {
    splicing_gff_converter_destroy(myconverter);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  return 0;
}

int splicing_iso_to_genomic_all(const splicing_gff_t *gff, size_t gene,
				int position, 
				const splicing_gff_converter_t *converter,
				splicing_vector_int_t *result) {

  size_t i;
  splicing_gff_converter_t vconverter, 
    *myconverter = (splicing_gff_converter_t*) converter;

  if (position < 1) { 
    SPLICING_ERROR("Invalid isoform coordinate, must the larger than zero", 
		   SPLICING_EINVAL);
  }

  if (!converter) { 
    myconverter=&vconverter;
    SPLICING_CHECK(splicing_gff_converter_init(gff, gene, myconverter));
    SPLICING_FINALLY(splicing_gff_converter_destroy, myconverter);
  }

  SPLICING_CHECK(splicing_vector_int_resize(result, myconverter->noiso));

  /* TODO: find impossible positions */
  for (i=0; i<myconverter->noiso; i++) {
    int ex;
    for (ex=VECTOR(myconverter->exidx)[i]; 
	 ex < VECTOR(myconverter->exidx)[i+1] && 
	   VECTOR(myconverter->exlim)[ex] <= position; 
	 ex++) ;
    if (ex < VECTOR(myconverter->exidx)[i+1]) {
      VECTOR(*result)[i] = position + VECTOR(myconverter->shift)[ex];
    } else {
      VECTOR(*result)[i] = -1;
    }
  }

  if (!converter) {
    splicing_gff_converter_destroy(myconverter);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  return 0;
}

int splicing_iso_to_genomic_1(const splicing_gff_t *gff, size_t gene,
			      int isoform, int position, 
			      const splicing_gff_converter_t *converter,
			      int *result) {

  int ex;
  splicing_gff_converter_t vconverter, 
    *myconverter = (splicing_gff_converter_t*) converter;

  if (position < 1) { 
    SPLICING_ERROR("Invalid isoform coordinate, must the larger than zero", 
		   SPLICING_EINVAL);
  }

  if (!converter) { 
    myconverter=&vconverter;
    SPLICING_CHECK(splicing_gff_converter_init(gff, gene, myconverter));
    SPLICING_FINALLY(splicing_gff_converter_destroy, myconverter);
  }

  /* TODO: find impossible positions */
  for (ex=VECTOR(myconverter->exidx)[isoform]; 
       ex < VECTOR(myconverter->exidx)[isoform+1] && 
	 VECTOR(myconverter->exlim)[ex] <= position; 
       ex++) ;
  if (ex < VECTOR(myconverter->exidx)[isoform+1]) {
    *result = position + VECTOR(myconverter->shift)[ex];
  } else {
    *result = -1;
  }

  if (!converter) {
    splicing_gff_converter_destroy(myconverter);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  return 0;
}

int splicing_genomic_to_iso(const splicing_gff_t *gff, size_t gene,
			    const splicing_vector_int_t *position, 
			    const splicing_gff_converter_t *converter,
			    splicing_matrix_int_t *isopos) {

  size_t r, i, noreads=splicing_vector_int_size(position);
  splicing_gff_converter_t vconverter, 
    *myconverter = (splicing_gff_converter_t*) converter;
  
  if (!converter) { 
    myconverter=&vconverter;
    SPLICING_CHECK(splicing_gff_converter_init(gff, gene, myconverter));
    SPLICING_FINALLY(splicing_gff_converter_destroy, myconverter);
  }

  SPLICING_CHECK(splicing_matrix_int_resize(isopos, myconverter->noiso, 
					    noreads));
  
  for (r=0; r<noreads; r++) {
    for (i=0; i<myconverter->noiso; i++) {
      size_t pos=VECTOR(*position)[r];
      size_t startpos=VECTOR(myconverter->exidx)[i];
      size_t endpos=VECTOR(myconverter->exidx)[i+1];
      int ex;
      for (ex=startpos; 
	   ex < endpos && VECTOR(myconverter->exend)[ex] < pos; 
	   ex++) ;
      if (ex < endpos && VECTOR(myconverter->exstart)[ex] <= pos && 
	  pos <= VECTOR(myconverter->exend)[ex]) {
	MATRIX(*isopos, i, r) = VECTOR(*position)[r] - 
	  VECTOR(myconverter->shift)[ex];
      } else { 
	MATRIX(*isopos, i, r) = -1;
      }
    }
  }

  if (!converter) { 
    splicing_gff_converter_destroy(myconverter);
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}

int splicing_genomic_to_iso_1(const splicing_gff_t *gff, size_t gene,
			      int isoform, int position, 
			      const splicing_gff_converter_t *converter,
			      int *result) {

  size_t startpos, endpos, ex;
  splicing_gff_converter_t vconverter, 
    *myconverter = (splicing_gff_converter_t*) converter;
  
  if (!converter) { 
    myconverter=&vconverter;
    SPLICING_CHECK(splicing_gff_converter_init(gff, gene, myconverter));
    SPLICING_FINALLY(splicing_gff_converter_destroy, myconverter);
  }

  startpos=VECTOR(myconverter->exidx)[isoform];
  endpos=VECTOR(myconverter->exidx)[isoform+1];
  for (ex=startpos; 
       ex < endpos && VECTOR(myconverter->exend)[ex] < position; 
       ex++) ;
  if (ex < endpos && VECTOR(myconverter->exstart)[ex] <= position && 
      position <= VECTOR(myconverter->exend)[ex]) {
    *result = position - VECTOR(myconverter->shift)[ex];
  } else { 
    *result = -1;
  }

  if (!converter) { 
    splicing_gff_converter_destroy(myconverter);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  return 0;
}

int splicing_genomic_to_iso_all(const splicing_gff_t *gff, size_t gene,
				int position, 
				const splicing_gff_converter_t *converter,
				splicing_vector_int_t *result) {

  int i;
  splicing_gff_converter_t vconverter, 
    *myconverter = (splicing_gff_converter_t*) converter;
  
  if (!converter) { 
    myconverter=&vconverter;
    SPLICING_CHECK(splicing_gff_converter_init(gff, gene, myconverter));
    SPLICING_FINALLY(splicing_gff_converter_destroy, myconverter);
  }

  SPLICING_CHECK(splicing_vector_int_resize(result, myconverter->noiso));
  
  for (i=0; i<myconverter->noiso; i++) {
    size_t startpos=VECTOR(myconverter->exidx)[i];
    size_t endpos=VECTOR(myconverter->exidx)[i+1];
    int ex;
    for (ex=startpos; 
	 ex < endpos && VECTOR(myconverter->exend)[ex] < position; 
	 ex++) ;
    if (ex < endpos && VECTOR(myconverter->exstart)[ex] <= position && 
	position <= VECTOR(myconverter->exend)[ex]) {
      VECTOR(*result)[i] = position - VECTOR(myconverter->shift)[ex];
    } else { 
      VECTOR(*result)[i] = -1;
    }
  }

  if (!converter) { 
    splicing_gff_converter_destroy(myconverter);
    SPLICING_FINALLY_CLEAN(1);
  }
  
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
