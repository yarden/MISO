
#include <R.h>
#include <Rdefines.h>

#include "splicing_random.h"
#include "splicing_error.h"
#include "Rsplicing.h"

void R_splicing_error_handler(const char *reason, const char *file, 
			      int line, int splicing_errno) {
  SPLICING_FINALLY_FREE();
  error("At %s:%i : %s, %s", file, line, reason,
	splicing_strerror(splicing_errno));
}

void R_splicing_warning_handler(const char *reason, const char *file,
				int line, int splicing_errno) {
  warning("At %s:%i : %s, %s", file, line, reason,
	  splicing_strerror(splicing_errno));
}

int R_splicing_begin() {
  splicing_set_error_handler(R_splicing_error_handler);
  splicing_set_warning_handler(R_splicing_warning_handler);
  GetRNGstate();
  return 0;
}

int R_splicing_end() {
  PutRNGstate();
  return 0;
}

SEXP R_splicing_finalizer() {
  SPLICING_FINALLY_FREE();
  return R_NilValue;
}

SEXP R_splicing_getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

SEXP R_splicing_vector_int_to_SEXP(const splicing_vector_int_t *v) {
  SEXP result;
  size_t n=splicing_vector_int_size(v);
  
  PROTECT(result=NEW_INTEGER(n));
  memcpy(INTEGER(result), VECTOR(*v), n*sizeof(int));
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_vector_to_SEXP(const splicing_vector_t *v) {
  SEXP result;
  size_t n=splicing_vector_size(v);
  
  PROTECT(result=NEW_NUMERIC(n));
  memcpy(REAL(result), VECTOR(*v), n*sizeof(double));
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_matrix_to_SEXP(const splicing_matrix_t *m) {
  int r = splicing_matrix_nrow(m);
  int c = splicing_matrix_ncol(m);
  int n = r * c;
  SEXP result, dim;
  
  PROTECT(result = NEW_NUMERIC(n));
  PROTECT(dim =  NEW_INTEGER(2));
  INTEGER(dim)[0] = r;
  INTEGER(dim)[1] = c;
  SET_DIM(result, dim);

  memcpy(REAL(result), &MATRIX(*m, 0, 0), n * sizeof(double));
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_matrix_int_to_SEXP(const splicing_matrix_int_t *m) {
  int r = splicing_matrix_int_nrow(m);
  int c = splicing_matrix_int_ncol(m);
  int n = r * c;
  SEXP result, dim;
  
  PROTECT(result = NEW_INTEGER(n));
  PROTECT(dim =  NEW_INTEGER(2));
  INTEGER(dim)[0] = r;
  INTEGER(dim)[1] = c;
  SET_DIM(result, dim);

  memcpy(INTEGER(result), &MATRIX(*m, 0, 0), n * sizeof(int));
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_vector_int_index_to_SEXP(const splicing_vector_int_t *v,
					 const splicing_vector_int_t *idx) {
  int i, n=splicing_vector_int_size(idx);
  int vn=splicing_vector_int_size(v);
  SEXP result;
  
  PROTECT(result=NEW_LIST(n));

  for (i=0; i<n; i++) {
    int j=VECTOR(*idx)[i];
    int x, k= i<n-1 ? VECTOR(*idx)[i+1] : vn;
    SEXP tmp;
    int *itmp;
    PROTECT(tmp=NEW_INTEGER(k-j));
    itmp=INTEGER(tmp);
    for (x=0; j<k; j++, x++) {
      itmp[x] = VECTOR(*v)[j];
    }
    SET_VECTOR_ELT(result, i, tmp);
    UNPROTECT(1);
  }
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_miso_rundata_to_SEXP(const splicing_miso_rundata_t *data) {
  SEXP result, names;
  
  PROTECT(result=NEW_LIST(6));
  SET_VECTOR_ELT(result, 0, ScalarInteger(data->noIso));
  SET_VECTOR_ELT(result, 1, ScalarInteger(data->noIters));
  SET_VECTOR_ELT(result, 2, ScalarInteger(data->noBurnIn));
  SET_VECTOR_ELT(result, 3, ScalarInteger(data->noLag));
  SET_VECTOR_ELT(result, 4, ScalarInteger(data->noAccepted));
  SET_VECTOR_ELT(result, 5, ScalarInteger(data->noRejected));
  
  PROTECT(names=NEW_STRING(6));
  SET_STRING_ELT(names, 0, mkChar("noIso"));
  SET_STRING_ELT(names, 1, mkChar("noIters"));
  SET_STRING_ELT(names, 2, mkChar("noBurnIn"));
  SET_STRING_ELT(names, 3, mkChar("noLag"));
  SET_STRING_ELT(names, 4, mkChar("noAccepted"));
  SET_STRING_ELT(names, 5, mkChar("noRejected"));
  
  SET_NAMES(result, names);
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_reads_to_SEXP(const splicing_reads_t *reads) {
  SEXP result, names, class;
  
  PROTECT(result=NEW_LIST(16));
  SET_VECTOR_ELT(result, 0, R_splicing_strvector_to_SEXP(&reads->chrname));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&reads->chrlen));
  SET_VECTOR_ELT(result, 2, R_splicing_vector_int_to_SEXP(&reads->chr));
  SET_VECTOR_ELT(result, 3, R_splicing_strvector_to_SEXP(&reads->qname));
  SET_VECTOR_ELT(result, 4, R_splicing_strvector_to_SEXP(&reads->cigar));
  SET_VECTOR_ELT(result, 5, R_splicing_vector_int_to_SEXP(&reads->position));
  SET_VECTOR_ELT(result, 6, R_splicing_vector_int_to_SEXP(&reads->flags));
  SET_VECTOR_ELT(result, 7, R_splicing_vector_int_to_SEXP(&reads->pairpos));
  SET_VECTOR_ELT(result, 8, ScalarInteger(reads->noPairs));
  SET_VECTOR_ELT(result, 9, ScalarInteger(reads->noSingles));
  SET_VECTOR_ELT(result, 10, ScalarLogical(reads->paired));
  SET_VECTOR_ELT(result, 11, R_splicing_vector_int_to_SEXP(&reads->mapq));
  SET_VECTOR_ELT(result, 12, R_splicing_vector_int_to_SEXP(&reads->rnext));
  SET_VECTOR_ELT(result, 13, R_splicing_vector_int_to_SEXP(&reads->tlen));
  SET_VECTOR_ELT(result, 14, R_splicing_strvector_to_SEXP(&reads->seq));
  SET_VECTOR_ELT(result, 15, R_splicing_strvector_to_SEXP(&reads->qual));

  PROTECT(names=NEW_CHARACTER(16));
  SET_STRING_ELT(names, 0, mkChar("chrname"));
  SET_STRING_ELT(names, 1, mkChar("chrlen"));
  SET_STRING_ELT(names, 2, mkChar("chr"));
  SET_STRING_ELT(names, 3, mkChar("qname"));
  SET_STRING_ELT(names, 4, mkChar("cigar"));
  SET_STRING_ELT(names, 5, mkChar("position"));
  SET_STRING_ELT(names, 6, mkChar("flag"));
  SET_STRING_ELT(names, 7, mkChar("pairpos"));
  SET_STRING_ELT(names, 8, mkChar("noPairs"));
  SET_STRING_ELT(names, 9, mkChar("noSingles"));
  SET_STRING_ELT(names, 10, mkChar("paired"));
  SET_STRING_ELT(names, 11, mkChar("mapq"));
  SET_STRING_ELT(names, 12, mkChar("rnext"));
  SET_STRING_ELT(names, 13, mkChar("tlen"));
  SET_STRING_ELT(names, 14, mkChar("seq"));
  SET_STRING_ELT(names, 15, mkChar("qual"));
  SET_NAMES(result, names);

  PROTECT(class=NEW_CHARACTER(1));
  SET_STRING_ELT(class, 0, mkChar("splicingSAM"));
  SET_CLASS(result, class);

  UNPROTECT(3);
  return result;
}

int R_splicing_SEXP_to_vector_int(SEXP pv, splicing_vector_int_t *v) {
  v->stor_begin = INTEGER(pv);
  v->stor_end = v->end = v->stor_begin + GET_LENGTH(pv);
  return 0;
}

int R_splicing_SEXP_to_vector(SEXP pv, splicing_vector_t *v) {
  v->stor_begin = REAL(pv);
  v->stor_end = v->end = v->stor_begin + GET_LENGTH(pv);
  return 0;
}

int R_splicing_SEXP_to_matrix(SEXP pm, splicing_matrix_t *m) {
  int *dim=INTEGER(GET_DIM(pm));
  R_splicing_SEXP_to_vector(pm, &m->data);
  m->nrow=dim[0];
  m->ncol=dim[1];
  return 0;
}

int R_splicing_SEXP_to_matrix_int(SEXP pm, splicing_matrix_int_t *m) {
  int *dim=INTEGER(GET_DIM(pm));
  R_splicing_SEXP_to_vector_int(pm, &m->data);
  m->nrow=dim[0];
  m->ncol=dim[1];
  return 0;
}

/* We construct a strvector by hand */

int R_splicing_SEXP_to_strvector(SEXP pv, splicing_strvector_t *v) {
  size_t i;
  v->size = v->asize = GET_LENGTH(pv);
  v->free = 0;
  v->table = (char **) R_alloc(v->size, sizeof(char*));
  for (i=0; i<v->size; i++) {
    v->table[i] = (char*) CHAR(STRING_ELT(pv, i));
  }

  return 1;
}

int R_splicing_SEXP_to_gff(SEXP pgff, splicing_gff_t *gff) {

  SEXP seqid, source, genes, transcripts, type, start, end, score, 
    strand, phase, ID, parent, seqid_str, source_str, attributes;

  seqid = R_splicing_getListElement(pgff, "seqid");
  source = R_splicing_getListElement(pgff, "source");
  genes = R_splicing_getListElement(pgff, "gid");
  transcripts = R_splicing_getListElement(pgff, "tid");
  type = R_splicing_getListElement(pgff, "type");
  start = R_splicing_getListElement(pgff, "start");
  end = R_splicing_getListElement(pgff, "end");
  score = R_splicing_getListElement(pgff, "score");
  strand = R_splicing_getListElement(pgff, "strand");
  phase = R_splicing_getListElement(pgff, "phase");
  ID = R_splicing_getListElement(pgff, "ID");
  parent = R_splicing_getListElement(pgff, "parent");
  seqid_str = R_splicing_getListElement(pgff, "seqid_str");
  source_str = R_splicing_getListElement(pgff, "source_str");
  attributes = R_splicing_getListElement(pgff, "attributes");

  R_splicing_SEXP_to_vector_int(seqid, &gff->seqid);
  R_splicing_SEXP_to_vector_int(source, &gff->source);
  R_splicing_SEXP_to_vector_int(genes, &gff->genes);
  R_splicing_SEXP_to_vector_int(transcripts, &gff->transcripts);
  R_splicing_SEXP_to_vector_int(type, &gff->type);
  R_splicing_SEXP_to_vector_int(start, &gff->start);
  R_splicing_SEXP_to_vector_int(end, &gff->end);
  R_splicing_SEXP_to_vector(score, &gff->score);
  R_splicing_SEXP_to_vector_int(strand, &gff->strand);
  R_splicing_SEXP_to_vector_int(phase, &gff->phase);
  R_splicing_SEXP_to_strvector(ID, &gff->ID);
  R_splicing_SEXP_to_vector_int(parent, &gff->parent);
  R_splicing_SEXP_to_strvector(seqid_str, &gff->seqids);
  R_splicing_SEXP_to_strvector(source_str, &gff->sources);

  gff->n = GET_LENGTH(start);

  return 0;
}

int R_splicing_SEXP_to_exons(SEXP pexons, splicing_vector_int_t *exons) {
  int i, p, noexons=GET_LENGTH(pexons);  
  int *tmp=(int*) R_alloc(noexons * 2, sizeof(int));
  splicing_vector_int_view(exons, tmp, noexons*2);

  for (i=0, p=0; i<noexons; i++) {
    int *v=INTEGER(VECTOR_ELT(pexons, i));
    VECTOR(*exons)[p++] = v[0];
    VECTOR(*exons)[p++] = v[1];
  }
  
  return 0;
}

int R_splicing_SEXP_to_isoforms(SEXP piso, splicing_vector_int_t *iso) {
  int i, p, veclen=0, noiso=GET_LENGTH(piso);
  int *tmp;
  
  for (i=0; i<noiso; i++) {
    veclen += GET_LENGTH(VECTOR_ELT(piso, i));
  }
  
  tmp = (int*) R_alloc(veclen+noiso, sizeof(int));
  splicing_vector_int_view(iso, tmp, veclen+noiso);
  
  for (i=0, p=0; i<noiso; i++) {
    int j, ilen=GET_LENGTH(VECTOR_ELT(piso, i));
    int *v=INTEGER(VECTOR_ELT(piso, i));
    for (j=0; j<ilen; j++) {
      VECTOR(*iso)[p++] = v[j]-1;
    }
    VECTOR(*iso)[p++] = -1;
  }
  
  return 0;
}

SEXP R_splicing_strvector_to_SEXP(const splicing_strvector_t *strv) {
  int i, n=splicing_strvector_size(strv);
  SEXP res;
  
  PROTECT(res=NEW_CHARACTER(n));
  for (i=0; i<n; i++) {
    SET_STRING_ELT(res, i, 
		   mkChar(splicing_strvector_get(strv, i)));
  }
  
  UNPROTECT(1);
  return res;
}

SEXP R_splicing_gff_to_SEXP(splicing_gff_t *gff) {
  SEXP result, class, names;

  PROTECT(result=NEW_LIST(15));
  SET_VECTOR_ELT(result, 0, R_splicing_strvector_to_SEXP(&gff->seqids));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&gff->seqid));
  SET_VECTOR_ELT(result, 2, R_splicing_strvector_to_SEXP(&gff->sources));
  SET_VECTOR_ELT(result, 3, R_splicing_vector_int_to_SEXP(&gff->source));
  SET_VECTOR_ELT(result, 4, R_splicing_vector_int_to_SEXP(&gff->type));
  SET_VECTOR_ELT(result, 5, R_splicing_vector_int_to_SEXP(&gff->start));
  SET_VECTOR_ELT(result, 6, R_splicing_vector_int_to_SEXP(&gff->end));
  SET_VECTOR_ELT(result, 7, R_splicing_vector_to_SEXP(&gff->score));
  SET_VECTOR_ELT(result, 8, R_splicing_vector_int_to_SEXP(&gff->strand));
  SET_VECTOR_ELT(result, 9, R_splicing_vector_int_to_SEXP(&gff->phase));
  SET_VECTOR_ELT(result, 10, R_NilValue); /* TODO: keep attributes */
  SET_VECTOR_ELT(result, 11, R_splicing_vector_int_to_SEXP(&gff->genes));
  SET_VECTOR_ELT(result, 12,
		 R_splicing_vector_int_to_SEXP(&gff->transcripts));
  SET_VECTOR_ELT(result, 13, R_splicing_strvector_to_SEXP(&gff->ID));
  SET_VECTOR_ELT(result, 14, R_splicing_vector_int_to_SEXP(&gff->parent));

  PROTECT(names=NEW_CHARACTER(15));
  SET_STRING_ELT(names, 0, mkChar("seqid_str"));
  SET_STRING_ELT(names, 1, mkChar("seqid"));
  SET_STRING_ELT(names, 2, mkChar("source_str"));  
  SET_STRING_ELT(names, 3, mkChar("source"));
  SET_STRING_ELT(names, 4, mkChar("type"));
  SET_STRING_ELT(names, 5, mkChar("start"));
  SET_STRING_ELT(names, 6, mkChar("end"));
  SET_STRING_ELT(names, 7, mkChar("score"));
  SET_STRING_ELT(names, 8, mkChar("strand"));
  SET_STRING_ELT(names, 9, mkChar("phase"));
  SET_STRING_ELT(names, 10, mkChar("attributes"));
  SET_STRING_ELT(names, 11, mkChar("gid"));
  SET_STRING_ELT(names, 12, mkChar("tid"));
  SET_STRING_ELT(names, 13, mkChar("ID"));
  SET_STRING_ELT(names, 14, mkChar("parent"));
  SET_NAMES(result, names);

  PROTECT(class=NEW_CHARACTER(1));
  SET_STRING_ELT(class, 0, mkChar("gff3"));
  SET_CLASS(result, class);

  UNPROTECT(3);
  return result;
}

SEXP R_splicing_read_gff(SEXP file) {
  const char *filename=CHAR(STRING_ELT(file, 0));
  FILE *input=fopen(filename, "r");
  splicing_gff_t gff;
  SEXP result;

  R_splicing_begin();
  
  splicing_gff_init(&gff, 50);
  splicing_gff_read(input, &gff);
  PROTECT(result = R_splicing_gff_to_SEXP(&gff));
  splicing_gff_destroy(&gff);
  fclose(input);
  
  R_splicing_end();

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_write_gff(SEXP pgff, SEXP file) {
  const char *filename=CHAR(STRING_ELT(file, 0));
  FILE *output=fopen(filename, "w");
  splicing_gff_t gff;
  SEXP result;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_gff_write(output, &gff);
  fclose(output);

  R_splicing_end();
  
  return R_NilValue;
}

SEXP R_splicing_matchIso(SEXP pgff, SEXP pgene, SEXP pposition, 
			 SEXP pcigar) {

  int i, noreads=GET_LENGTH(pcigar);
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t exstart, exend, position;
  splicing_gff_t gff;
  splicing_matrix_t res;
  SEXP result;
  const char **cigarstr;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  splicing_matrix_init(&res, 0, 0);
  cigarstr = (const char**) R_alloc(noreads, sizeof(char*));
  
  for (i=0; i<noreads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }

  splicing_matchIso(&gff, gene, &position, cigarstr, &res);

  PROTECT(result = R_splicing_matrix_to_SEXP(&res));

  splicing_matrix_destroy(&res);

  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_parse_cigar(SEXP cigar) {
  size_t i, n=GET_LENGTH(cigar);
  SEXP result, names; 
  splicing_vector_int_t numcigar, cigaridx;
  const char **strcigar;
  
  R_splicing_begin();
  
  PROTECT(result = NEW_LIST(2));
  PROTECT(names = NEW_CHARACTER(2));
  SET_STRING_ELT(names, 0, mkChar("cigar"));
  SET_STRING_ELT(names, 1, mkChar("cigaridx"));
  SET_NAMES(result, names);
  
  splicing_vector_int_init(&numcigar, 0);
  splicing_vector_int_init(&cigaridx, 0);
  strcigar = (const char**) R_alloc(n, sizeof(char*));

  for (i=0; i<n; i++) {
    strcigar[i] = CHAR(STRING_ELT(cigar, i));
  }
  
  splicing_parse_cigar(strcigar, n, &numcigar, &cigaridx);
  
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&numcigar));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&cigaridx));
  
  splicing_vector_int_destroy(&numcigar);
  splicing_vector_int_destroy(&cigaridx);  
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_order_matches(SEXP pmatches) {
  SEXP result; 
  splicing_vector_int_t order;
  splicing_matrix_t matches;

  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix(pmatches, &matches);
  splicing_vector_int_init(&order, 0); 
  
  splicing_order_matches(&matches, &order);
  
  PROTECT(result=R_splicing_vector_int_to_SEXP(&order));
  
  splicing_vector_int_destroy(&order);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_miso(SEXP pgff, SEXP pgene, SEXP preads, SEXP preadLength, 
		     SEXP pnoIterations, SEXP pnoBurnIn, SEXP pnoLag,
		     SEXP phyperp) {
  
  size_t gene=INTEGER(pgene)[0]-1;
  SEXP result, names, class;
  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_gff_t gff;
  splicing_vector_int_t position;
  const char **cigarstr;
  SEXP cigar=R_splicing_getListElement(preads, "cigar");
  int readLength=INTEGER(preadLength)[0];
  int noIterations=INTEGER(pnoIterations)[0];
  int noBurnIn=INTEGER(pnoBurnIn)[0];
  int noLag=INTEGER(pnoLag)[0];
  splicing_vector_t hyperp;
  int i, noReads=GET_LENGTH(cigar);
  splicing_matrix_t match_matrix;
  splicing_matrix_t class_templates;
  splicing_vector_t class_counts;
  splicing_miso_rundata_t rundata;

  R_splicing_begin();
  
  splicing_matrix_init(&samples, 0, 0);
  splicing_vector_init(&logLik, 0);
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&class_templates, 0, 0);
  splicing_vector_init(&class_counts, 0);

  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(R_splicing_getListElement(preads, "position"),
				&position);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  R_splicing_SEXP_to_vector(phyperp, &hyperp);

  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(cigar, i));
  }

  splicing_miso(&gff, gene, &position, cigarstr, readLength, noIterations,
		noBurnIn, noLag, &hyperp, &samples, &logLik, 
		&match_matrix, &class_templates, &class_counts, 
		/*assignment=*/ 0, &rundata);
  
  PROTECT(result=NEW_LIST(6));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&samples));
  splicing_matrix_destroy(&samples);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&logLik));
  splicing_vector_destroy(&logLik);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 3, R_splicing_matrix_to_SEXP(&class_templates));
  splicing_matrix_destroy(&class_templates);
  SET_VECTOR_ELT(result, 4, R_splicing_vector_to_SEXP(&class_counts));
  splicing_vector_destroy(&class_counts);
  SET_VECTOR_ELT(result, 5, R_splicing_miso_rundata_to_SEXP(&rundata));

  PROTECT(names=NEW_CHARACTER(6));
  SET_STRING_ELT(names, 0, mkChar("samples"));
  SET_STRING_ELT(names, 1, mkChar("logLik"));
  SET_STRING_ELT(names, 2, mkChar("matchMatrix"));
  SET_STRING_ELT(names, 3, mkChar("classTemplates"));
  SET_STRING_ELT(names, 4, mkChar("classCounts"));
  SET_STRING_ELT(names, 5, mkChar("runData"));
  SET_NAMES(result, names);

  PROTECT(class=ScalarString(mkChar("MISOresult")));
  SET_CLASS(result, class);

  R_splicing_end();
  
  UNPROTECT(3);
  return result;
}

SEXP R_splicing_miso_trinity(SEXP pmatch_matrix, SEXP pisolen, 
			     SEXP preadlength, SEXP pnoiterations,
			     SEXP pnoburnin, SEXP pnolag, SEXP phyperp) {
  
  SEXP result, names, class;
  splicing_matrix_t match_matrix;
  splicing_vector_int_t isolen;
  int readlength=INTEGER(preadlength)[0];
  int noiterations=INTEGER(pnoiterations)[0];
  int noburnin=INTEGER(pnoburnin)[0];
  int nolag=INTEGER(pnolag)[0];
  splicing_vector_t hyperp;
  
  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_matrix_t class_templates;
  splicing_vector_t class_counts;
  splicing_miso_rundata_t rundata;

  R_splicing_begin();
  
  splicing_matrix_init(&samples, 0, 0);
  splicing_vector_init(&logLik, 0);
  splicing_matrix_init(&class_templates, 0, 0);
  splicing_vector_init(&class_counts, 0);

  R_splicing_SEXP_to_matrix(pmatch_matrix, &match_matrix);
  R_splicing_SEXP_to_vector_int(pisolen, &isolen);
  R_splicing_SEXP_to_vector(phyperp, &hyperp);
  
  splicing_miso_trinity(&match_matrix, &isolen, readlength, noiterations,
			noburnin, nolag, &hyperp, &samples, &logLik, 
			&class_templates, &class_counts, /*assignment=*/ 0,
			&rundata);
  
  PROTECT(result=NEW_LIST(5));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&samples));
  splicing_matrix_destroy(&samples);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&logLik));
  splicing_vector_destroy(&logLik);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&class_templates));
  splicing_matrix_destroy(&class_templates);
  SET_VECTOR_ELT(result, 3, R_splicing_vector_to_SEXP(&class_counts));
  splicing_vector_destroy(&class_counts);
  SET_VECTOR_ELT(result, 4, R_splicing_miso_rundata_to_SEXP(&rundata));

  PROTECT(names=NEW_CHARACTER(5));
  SET_STRING_ELT(names, 0, mkChar("samples"));
  SET_STRING_ELT(names, 1, mkChar("logLik"));
  SET_STRING_ELT(names, 2, mkChar("classTemplates"));
  SET_STRING_ELT(names, 3, mkChar("classCounts"));
  SET_STRING_ELT(names, 4, mkChar("runData"));
  SET_NAMES(result, names);

  PROTECT(class=ScalarString(mkChar("MISOresult")));
  SET_CLASS(result, class);

  R_splicing_end();
  
  UNPROTECT(3);
  return result;
}  

SEXP R_splicing_miso_paired(SEXP pgff, SEXP pgene, SEXP preads, 
			    SEXP preadLength, SEXP pnoIterations, 
			    SEXP pnoBurnIn, SEXP pnoLag, SEXP phyperp, 
			    SEXP pfragmentProb, SEXP pfragmentStart, 
			    SEXP pnormalMean, SEXP pnormalVar, 
			    SEXP pnumDevs) {
  
  size_t gene=INTEGER(pgene)[0]-1;
  SEXP result, names, class;
  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_gff_t gff;
  splicing_vector_int_t position;
  const char **cigarstr;
  SEXP cigar=R_splicing_getListElement(preads, "cigar");
  int readLength=INTEGER(preadLength)[0];
  int noIterations=INTEGER(pnoIterations)[0];
  int noBurnIn=INTEGER(pnoBurnIn)[0];
  int noLag=INTEGER(pnoLag)[0];
  splicing_vector_t hyperp;
  splicing_vector_t fragmentProb;
  int fragmentStart=INTEGER(pfragmentStart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  int i, noReads=GET_LENGTH(cigar);
  splicing_matrix_t match_matrix;
  splicing_matrix_t class_templates;
  splicing_vector_t class_counts;
  splicing_miso_rundata_t rundata;

  R_splicing_begin();
  
  splicing_matrix_init(&samples, 0, 0);
  splicing_vector_init(&logLik, 0);
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&class_templates, 0, 0);
  splicing_vector_init(&class_counts, 0);

  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(R_splicing_getListElement(preads, "position"),
				&position);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  R_splicing_SEXP_to_vector(phyperp, &hyperp);

  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(cigar, i));
  }
  
  if (!isNull(pfragmentProb)) {
    R_splicing_SEXP_to_vector(pfragmentProb, &fragmentProb);
  }

  splicing_miso_paired(&gff, gene, &position, cigarstr, readLength,
		       noIterations, noBurnIn, noLag, &hyperp, 
		       isNull(pfragmentProb) ? 0 : &fragmentProb, 
		       fragmentStart, normalMean, normalVar, numDevs,
		       &samples, &logLik, &match_matrix, &class_templates,
		       &class_counts, /*assignment=*/ 0, &rundata);
  
  PROTECT(result=NEW_LIST(6));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&samples));
  splicing_matrix_destroy(&samples);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&logLik));
  splicing_vector_destroy(&logLik);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 3, R_splicing_matrix_to_SEXP(&class_templates));
  splicing_matrix_destroy(&class_templates);
  SET_VECTOR_ELT(result, 4, R_splicing_vector_to_SEXP(&class_counts));
  splicing_vector_destroy(&class_counts);
  SET_VECTOR_ELT(result, 5, R_splicing_miso_rundata_to_SEXP(&rundata));

  PROTECT(names=NEW_CHARACTER(6));
  SET_STRING_ELT(names, 0, mkChar("samples"));
  SET_STRING_ELT(names, 1, mkChar("logLik"));
  SET_STRING_ELT(names, 2, mkChar("matchMatrix"));
  SET_STRING_ELT(names, 3, mkChar("classTemplates"));
  SET_STRING_ELT(names, 4, mkChar("classCounts"));
  SET_STRING_ELT(names, 5, mkChar("runData"));
  SET_NAMES(result, names);

  PROTECT(class=ScalarString(mkChar("MISOresult")));
  SET_CLASS(result, class);

  R_splicing_end();
  
  UNPROTECT(3);
  return result;
}

SEXP R_splicing_miso_paired_trinity(SEXP pmatch_matrix, 
				    SEXP pfragment_length,
				    SEXP pisolen,
				    SEXP preadLength, SEXP pnoIterations, 
				    SEXP pnoBurnIn, SEXP pnoLag, 
				    SEXP phyperp, 
				    SEXP pfragmentProb, SEXP pfragmentStart, 
				    SEXP pnormalMean, SEXP pnormalVar, 
				    SEXP pnumDevs) {
  
  SEXP result, names, class;

  splicing_matrix_t match_matrix;
  splicing_matrix_int_t fragment_length;
  splicing_vector_int_t isolen;
  int readLength=INTEGER(preadLength)[0];
  int noIterations=INTEGER(pnoIterations)[0];
  int noBurnIn=INTEGER(pnoBurnIn)[0];
  int noLag=INTEGER(pnoLag)[0];
  splicing_vector_t hyperp;
  splicing_vector_t fragmentProb;
  int fragmentStart=INTEGER(pfragmentStart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];

  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_matrix_t class_templates;
  splicing_vector_t class_counts;
  splicing_miso_rundata_t rundata;

  R_splicing_begin();
  
  splicing_matrix_init(&samples, 0, 0);
  splicing_vector_init(&logLik, 0);
  splicing_matrix_init(&class_templates, 0, 0);
  splicing_vector_init(&class_counts, 0);

  R_splicing_SEXP_to_matrix(pmatch_matrix, &match_matrix);
  R_splicing_SEXP_to_matrix_int(pfragment_length, &fragment_length);
  R_splicing_SEXP_to_vector_int(pisolen, &isolen);
  R_splicing_SEXP_to_vector(phyperp, &hyperp);

  if (!isNull(pfragmentProb)) {
    R_splicing_SEXP_to_vector(pfragmentProb, &fragmentProb);
  }

  splicing_miso_paired_trinity(&match_matrix, &fragment_length, &isolen, 
			       readLength, noIterations, noBurnIn, noLag, 
			       &hyperp, 
			       isNull(pfragmentProb) ? 0 : &fragmentProb, 
			       fragmentStart, normalMean, normalVar, numDevs,
			       &samples, &logLik, &class_templates,
			       &class_counts, /*assignment=*/ 0, &rundata);
  
  PROTECT(result=NEW_LIST(5));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&samples));
  splicing_matrix_destroy(&samples);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&logLik));
  splicing_vector_destroy(&logLik);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&class_templates));
  splicing_matrix_destroy(&class_templates);
  SET_VECTOR_ELT(result, 3, R_splicing_vector_to_SEXP(&class_counts));
  splicing_vector_destroy(&class_counts);
  SET_VECTOR_ELT(result, 4, R_splicing_miso_rundata_to_SEXP(&rundata));

  PROTECT(names=NEW_CHARACTER(5));
  SET_STRING_ELT(names, 0, mkChar("samples"));
  SET_STRING_ELT(names, 1, mkChar("logLik"));
  SET_STRING_ELT(names, 2, mkChar("classTemplates"));
  SET_STRING_ELT(names, 3, mkChar("classCounts"));
  SET_STRING_ELT(names, 4, mkChar("runData"));
  SET_NAMES(result, names);

  PROTECT(class=ScalarString(mkChar("MISOresult")));
  SET_CLASS(result, class);

  R_splicing_end();
  
  UNPROTECT(3);
  return result;
}

SEXP R_splicing_reassign_samples(SEXP pmatches, SEXP pmatch_order, 
				 SEXP ppsi, SEXP pnoiso) {
  SEXP result; 
  splicing_matrix_t matches;
  splicing_vector_int_t match_order;
  splicing_vector_t psi;
  int noiso=INTEGER(pnoiso)[0];
  splicing_vector_int_t cresult;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix(pmatches, &matches);
  R_splicing_SEXP_to_vector_int(pmatch_order, &match_order);
  R_splicing_SEXP_to_vector(ppsi, &psi);
  splicing_vector_int_init(&cresult, 0);

  splicing_reassign_samples(&matches, &match_order, &psi, noiso, &cresult);

  PROTECT(result=R_splicing_vector_int_to_SEXP(&cresult));
  
  splicing_vector_int_destroy(&cresult);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_mvplogisnorm(SEXP ptheta, SEXP pmu, SEXP psigma, 
			     SEXP plen) {
  SEXP result;
  splicing_vector_t theta, mu;
  double sigma=REAL(psigma)[0];
  int len=INTEGER(plen)[0];
  double score;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(ptheta, &theta);
  R_splicing_SEXP_to_vector(pmu, &mu);
  
  splicing_mvplogisnorm(&theta, &mu, sigma, len, &score);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = sigma;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_score_iso(SEXP ppsi, SEXP pnoiso, SEXP passignment, 
			  SEXP pnoreads, SEXP peffisolen) {
  SEXP result;
  splicing_vector_t psi;
  int noiso=INTEGER(pnoiso)[0];
  splicing_vector_int_t assignment;
  int noreads=INTEGER(pnoreads)[0];
  splicing_vector_int_t effisolen;
  double res;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(ppsi, &psi);
  R_splicing_SEXP_to_vector_int(passignment, &assignment);
  R_splicing_SEXP_to_vector_int(peffisolen, &effisolen);
  
  splicing_score_iso(&psi, noiso, &assignment, noreads, &effisolen, &res);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = res;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_ldirichlet(SEXP px, SEXP palpha, SEXP plen) {
  SEXP result;
  splicing_vector_t x, alpha;
  int len=INTEGER(plen)[0];
  double res;

  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(px, &x);
  R_splicing_SEXP_to_vector(palpha, &alpha);
 
  splicing_ldirichlet(&x, &alpha, len, &res);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = res;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_mvrnorm(SEXP pmu, SEXP psigma, SEXP plen) {
  SEXP result;
  splicing_vector_t mu, resalpha;
  double sigma=REAL(psigma)[0];
  int len=INTEGER(plen)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(pmu, &mu);
  splicing_vector_init(&resalpha, 0);

  splicing_mvrnorm(&mu, sigma, &resalpha, len);

  PROTECT(result=R_splicing_vector_to_SEXP(&resalpha));
  
  splicing_vector_destroy(&resalpha);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_logit_inv(SEXP px, SEXP plen) {
  SEXP result;
  splicing_vector_t x;
  splicing_vector_t res;
  int len=INTEGER(plen)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(px, &x);
  splicing_vector_init(&res, 0);
  
  splicing_logit_inv(&x, &res, len);
  
  PROTECT(result=R_splicing_vector_to_SEXP(&res));
  
  splicing_vector_destroy(&res);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_score_joint(SEXP passignment, SEXP pnoreads, SEXP ppsi,
			    SEXP phyper, SEXP peffisolen, SEXP pisoscores) {
  
  SEXP result;
  splicing_vector_int_t assignment, effisolen;
  splicing_vector_t psi, hyper, isoscores;
  int noreads=INTEGER(pnoreads)[0];
  double score;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector_int(passignment, &assignment);
  R_splicing_SEXP_to_vector_int(peffisolen, &effisolen);
  R_splicing_SEXP_to_vector(ppsi, &psi);
  R_splicing_SEXP_to_vector(phyper, &hyper);
  
  splicing_score_joint(&assignment, noreads, &psi, &hyper, &effisolen,
		       &isoscores, &score);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = score;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_drift_proposal(SEXP pmode, SEXP ppsi, SEXP palpha, 
			       SEXP psigma, SEXP potherpsi, SEXP potheralpha,
			       SEXP pnoiso) {
  int mode=INTEGER(pmode)[0];
  splicing_vector_t psi, alpha, otherpsi, otheralpha, respsi, resalpha;
  double sigma=REAL(psigma)[0];
  double ressigma, resscore;
  int noiso=INTEGER(pnoiso)[0];
  SEXP result, names;

  R_splicing_begin();
  
  switch (mode) {

  case 0:
    splicing_vector_init(&respsi, 0);
    splicing_vector_init(&resalpha, 0);
    splicing_drift_proposal(0, 0, 0, 0, 0, 0, noiso, &respsi, &resalpha, 
			    &ressigma, 0);
    PROTECT(result=NEW_LIST(3));
    SET_VECTOR_ELT(result, 0, R_splicing_vector_to_SEXP(&respsi));
    SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&resalpha));
    SET_VECTOR_ELT(result, 2, ScalarReal(ressigma));
    PROTECT(names=NEW_CHARACTER(3));
    SET_STRING_ELT(names, 0, mkChar("psi"));
    SET_STRING_ELT(names, 1, mkChar("alpha"));
    SET_STRING_ELT(names, 2, mkChar("sigma"));
    SET_NAMES(result, names);
    splicing_vector_destroy(&respsi);
    splicing_vector_destroy(&resalpha);
    UNPROTECT(2);
    break;

  case 1:
    R_splicing_SEXP_to_vector(ppsi, &psi);
    R_splicing_SEXP_to_vector(palpha, &alpha);
    splicing_vector_init(&respsi, 0);
    splicing_vector_init(&resalpha, 0);
    splicing_drift_proposal(1, &psi, &alpha, sigma, 0, 0, noiso, 
			    &respsi, &resalpha, 0, 0);
    PROTECT(result=NEW_LIST(2));
    SET_VECTOR_ELT(result, 0, R_splicing_vector_to_SEXP(&respsi));
    SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&resalpha));
    PROTECT(names=NEW_CHARACTER(2));
    SET_STRING_ELT(names, 0, mkChar("psi"));
    SET_STRING_ELT(names, 1, mkChar("alpha"));
    SET_NAMES(result, names);
    splicing_vector_destroy(&respsi);
    splicing_vector_destroy(&resalpha);
    UNPROTECT(2);
    break;

  case 2:
    R_splicing_SEXP_to_vector(ppsi, &psi);
    R_splicing_SEXP_to_vector(palpha, &alpha);
    splicing_drift_proposal(2, &psi, &alpha, sigma, &otherpsi, &otheralpha,
			    noiso, 0, 0, 0, &resscore);
    PROTECT(result=ScalarReal(resscore));
    UNPROTECT(1);
    break;
  }

  R_splicing_end();
  
  return result;
}

SEXP R_splicing_metropolis_hastings_ratio(SEXP passignment, SEXP pnoreads,
					  SEXP ppsiNew, SEXP palphaNew,
					  SEXP ppsi, SEXP palpha, 
					  SEXP psigma, SEXP pnoiso, 
					  SEXP peffisolen, SEXP phyperp,
					  SEXP pisoscores, SEXP pfull) {

  SEXP result, names; 
  splicing_vector_int_t assignment, effisolen;
  splicing_vector_t psiNew, alphaNew, psi, alpha, hyperp, isoscores;
  int noreads=INTEGER(pnoreads)[0]; 
  int noiso=INTEGER(pnoiso)[0];
  double sigma=REAL(psigma)[0];
  int full=INTEGER(pfull)[0];
  double acceptP, pcJS, ppJS;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector_int(passignment, &assignment);
  R_splicing_SEXP_to_vector(ppsiNew, &psiNew);
  R_splicing_SEXP_to_vector(palphaNew, &alphaNew);
  R_splicing_SEXP_to_vector(ppsi, &psi);
  R_splicing_SEXP_to_vector(palpha, &alpha);
  R_splicing_SEXP_to_vector_int(peffisolen, &effisolen);
  R_splicing_SEXP_to_vector(phyperp, &hyperp);
  R_splicing_SEXP_to_vector(pisoscores, &isoscores);
  
  splicing_metropolis_hastings_ratio(&assignment, noreads, &psiNew, 
				     &alphaNew, &psi, &alpha, sigma, noiso,
				     &effisolen, &hyperp, &isoscores, full,
				     &acceptP, &pcJS, &ppJS);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, ScalarReal(acceptP));
  SET_VECTOR_ELT(result, 1, ScalarReal(pcJS));
  SET_VECTOR_ELT(result, 2, ScalarReal(ppJS));
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(result, 0, mkChar("acceptP"));
  SET_STRING_ELT(result, 1, mkChar("pcJS"));
  SET_STRING_ELT(result, 2, mkChar("ppJS"));
  SET_NAMES(result, names);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;  
}

SEXP R_splicing_gff_exon_start_end(SEXP pgff, SEXP pgene) {
  SEXP result, names;
  int gene=INTEGER(pgene)[0]-1;
  splicing_gff_t gff;
  splicing_vector_int_t start, end, idx;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&start, 0);
  splicing_vector_int_init(&end, 0);
  splicing_vector_int_init(&idx, 0);
  
  splicing_gff_exon_start_end(&gff, &start, &end, &idx, gene);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&start));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&end));
  SET_VECTOR_ELT(result, 2, R_splicing_vector_int_to_SEXP(&idx));
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("start"));
  SET_STRING_ELT(names, 1, mkChar("end"));
  SET_STRING_ELT(names, 2, mkChar("idx"));
  SET_NAMES(result, names);

  splicing_vector_int_destroy(&idx);
  splicing_vector_int_destroy(&end);
  splicing_vector_int_destroy(&start);

  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}
  
					  
SEXP R_splicing_numeric_cigar(SEXP pexstart, SEXP pexend, SEXP pexidx, 
			      SEXP pnoiso, SEXP pgenestart) {
  SEXP result;
  splicing_vector_int_t exstart, exend, exidx, res;
  int noiso=INTEGER(pnoiso)[0];
  int genestart=INTEGER(pgenestart)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector_int(pexstart, &exstart);
  R_splicing_SEXP_to_vector_int(pexend, &exend);  
  R_splicing_SEXP_to_vector_int(pexidx, &exidx);  
  splicing_vector_int_init(&res, 0);

  splicing_numeric_cigar(&exstart, &exend, &exidx, noiso, genestart, &res);
  
  PROTECT(result=R_splicing_vector_int_to_SEXP(&res));
  
  splicing_vector_int_destroy(&res);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_assignment_matrix(SEXP pgff, SEXP pgene, SEXP preadlength) {
  SEXP result;
  splicing_matrix_t matrix;
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readlength=INTEGER(preadlength)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_matrix_init(&matrix, 0, 0);
  
  splicing_assignment_matrix(&gff, gene, readlength, &matrix);
  
  PROTECT(result=R_splicing_matrix_to_SEXP(&matrix));
  
  splicing_matrix_destroy(&matrix);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_create_gene(SEXP pexons, SEXP pisoforms, SEXP pid, 
			    SEXP pseqid, SEXP psource, SEXP pstrand) {
  SEXP result;
  splicing_vector_int_t exons, isoforms;
  const char *id=CHAR(STRING_ELT(pid, 0));
  const char *seqid=CHAR(STRING_ELT(pseqid, 0));
  const char *source=CHAR(STRING_ELT(psource, 0));
  int strand=INTEGER(pstrand)[0];
  splicing_gff_t gff;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_exons(pexons, &exons);
  R_splicing_SEXP_to_isoforms(pisoforms, &isoforms);
  splicing_gff_init(&gff, 0);
  
  splicing_create_gene(&exons, &isoforms, id, seqid, source, strand, &gff);
  
  PROTECT(result=R_splicing_gff_to_SEXP(&gff));
  
  splicing_gff_destroy(&gff);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_simulate_reads(SEXP pgff, SEXP pgene, SEXP pexpression,
			       SEXP pnoreads, SEXP preadLength) {
  SEXP result, names;
  size_t gene=INTEGER(pgene)[0]-1;
  int noreads=INTEGER(pnoreads)[0];
  int readlength=INTEGER(preadLength)[0];
  splicing_gff_t gff;
  splicing_vector_t expression;
  splicing_vector_int_t isoform, position;
  splicing_strvector_t cigar;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector(pexpression, &expression);
  splicing_vector_int_init(&isoform, 0);
  splicing_vector_int_init(&position, 0);
  splicing_strvector_init(&cigar, 0);

  splicing_simulate_reads(&gff, gene, &expression, noreads, readlength,
			  &isoform, &position, &cigar);

  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&isoform));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&position));
  SET_VECTOR_ELT(result, 2, R_splicing_strvector_to_SEXP(&cigar));
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("isoform"));
  SET_STRING_ELT(names, 1, mkChar("position"));
  SET_STRING_ELT(names, 2, mkChar("cigar"));
  SET_NAMES(result, names);
  
  splicing_strvector_destroy(&cigar);
  splicing_vector_int_destroy(&position);
  splicing_vector_int_destroy(&isoform);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_simulate_paired_reads(SEXP pgff, SEXP pgene, SEXP pexpression,
				      SEXP pnoReads, SEXP preadLength,
				      SEXP pfragmentProb, SEXP pfragmentStart, 
				      SEXP pnormalMean, SEXP pnormalVar,
				      SEXP pnumDevs) {
  
  SEXP result, names;
  size_t gene=INTEGER(pgene)[0]-1;
  int noreads=INTEGER(pnoReads)[0];
  int readLength=INTEGER(preadLength)[0];
  splicing_gff_t gff;
  splicing_vector_t expression, fragmentProb;
  int fragmentStart=INTEGER(pfragmentStart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  splicing_vector_int_t isoform, position;
  splicing_strvector_t cigar;
  splicing_vector_t sampleprob;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector(pexpression, &expression);
  if (!isNull(pfragmentProb)) { 
    R_splicing_SEXP_to_vector(pfragmentProb, &fragmentProb);
  }
  splicing_vector_int_init(&isoform, 0);
  splicing_vector_int_init(&position, 0);
  splicing_strvector_init(&cigar, 0);
  splicing_vector_init(&sampleprob, 0);
  
  splicing_simulate_paired_reads(&gff, gene, &expression, noreads, 
				 readLength, 
				 isNull(pfragmentProb) ? 0 : &fragmentProb, 
				 fragmentStart, normalMean, normalVar, 
				 numDevs, &isoform, &position, &cigar, 
				 &sampleprob);
  
  PROTECT(result=NEW_LIST(4));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&isoform));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&position));
  SET_VECTOR_ELT(result, 2, R_splicing_strvector_to_SEXP(&cigar));
  SET_VECTOR_ELT(result, 3, R_splicing_vector_to_SEXP(&sampleprob));
  PROTECT(names=NEW_CHARACTER(4));
  SET_STRING_ELT(names, 0, mkChar("isoform"));
  SET_STRING_ELT(names, 1, mkChar("position"));
  SET_STRING_ELT(names, 2, mkChar("cigar"));
  SET_STRING_ELT(names, 3, mkChar("sampleProb"));
  SET_NAMES(result, names);
  
  splicing_strvector_destroy(&cigar);
  splicing_vector_int_destroy(&position);
  splicing_vector_int_destroy(&isoform);
  splicing_vector_destroy(&sampleprob);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_gff_noiso(SEXP pgff) {
  SEXP result;
  splicing_gff_t gff;
  splicing_vector_int_t noiso;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&noiso, 0);
  
  splicing_gff_noiso(&gff, &noiso);
  
  PROTECT(result=R_splicing_vector_int_to_SEXP(&noiso));
  
  splicing_vector_int_destroy(&noiso);
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_gff_isolength(SEXP pgff) {
  SEXP result;
  splicing_gff_t gff;
  splicing_vector_int_t isolength, isolength_idx;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&isolength, 0);
  splicing_vector_int_init(&isolength_idx, 0);
  
  splicing_gff_isolength(&gff, &isolength, &isolength_idx);
  
  PROTECT(result=R_splicing_vector_int_index_to_SEXP(&isolength, 
						     &isolength_idx));
  
  splicing_vector_int_destroy(&isolength);
  splicing_vector_int_destroy(&isolength_idx);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_solve_gene(SEXP pgff, SEXP pgene, SEXP preadLength, 
			   SEXP pposition, SEXP pcigar) {

  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readLength=INTEGER(preadLength)[0];
  splicing_vector_int_t position;
  const char **cigarstr;
  int i, noReads=GET_LENGTH(pcigar);
  splicing_matrix_t match_matrix, assignment_matrix;
  splicing_vector_t expression;
  double rnorm;
  SEXP result, names;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&assignment_matrix, 0, 0);
  splicing_vector_init(&expression, 0);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }
  
  splicing_solve_gene(&gff, gene, readLength, &position, cigarstr, 
		      &match_matrix, &assignment_matrix, &expression);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 1, R_splicing_matrix_to_SEXP(&assignment_matrix));
  splicing_matrix_destroy(&assignment_matrix);
  SET_VECTOR_ELT(result, 2, R_splicing_vector_to_SEXP(&expression));
  splicing_vector_destroy(&expression);
  
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("match"));
  SET_STRING_ELT(names, 1, mkChar("assignment"));
  SET_STRING_ELT(names, 2, mkChar("expression"));
  SET_NAMES(result, names); 
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;		 
}

SEXP R_splicing_solve_gene_paired(SEXP pgff, SEXP pgene, SEXP preadLength, 
				  SEXP pposition, SEXP pcigar, 
				  SEXP pfragmentprob, SEXP pfragmentstart, 
				  SEXP pnormalMean, SEXP pnormalVar, 
				  SEXP pnumDevs) {

  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readLength=INTEGER(preadLength)[0];
  splicing_vector_int_t position;
  const char **cigarstr;
  int i, noReads=GET_LENGTH(pcigar);
  splicing_vector_t fragmentprob;
  int fragmentstart=INTEGER(pfragmentstart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  splicing_matrix_t match_matrix, assignment_matrix;
  splicing_vector_t expression;
  double rnorm;
  SEXP result, names;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  if (!isNull(pfragmentprob)) { 
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&assignment_matrix, 0, 0);
  splicing_vector_init(&expression, 0);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }
  
  splicing_solve_gene_paired(&gff, gene, readLength, &position, cigarstr, 
			     isNull(pfragmentprob) ? 0 : &fragmentprob, 
			     fragmentstart, normalMean, normalVar, numDevs,
			     &match_matrix, &assignment_matrix, &expression);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 1, R_splicing_matrix_to_SEXP(&assignment_matrix));
  splicing_matrix_destroy(&assignment_matrix);
  SET_VECTOR_ELT(result, 2, R_splicing_vector_to_SEXP(&expression));
  splicing_vector_destroy(&expression);
  
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("match"));
  SET_STRING_ELT(names, 1, mkChar("assignment"));
  SET_STRING_ELT(names, 2, mkChar("expression"));
  SET_NAMES(result, names); 
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;		 
}

SEXP R_splicing_iso_to_genomic(SEXP pgff, SEXP pgene, SEXP pisoform,
			       SEXP pposition) {
  
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t isoform, position, newposition;
  SEXP result;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pisoform, &isoform);

  PROTECT(result=NEW_INTEGER(GET_LENGTH(pposition)));
  memcpy(INTEGER(result), INTEGER(pposition), 
	 sizeof(int) * GET_LENGTH(pposition));
  splicing_vector_int_view(&newposition, INTEGER(result), 
			   GET_LENGTH(pposition));
  
  splicing_iso_to_genomic(&gff, gene, &isoform, 0, 0, 0, &newposition);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_genomic_to_iso(SEXP pgff, SEXP pgene, SEXP pposition) {

  SEXP result;
  
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t position;
  splicing_matrix_int_t isopos;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  splicing_matrix_int_init(&isopos, 0, 0);

  splicing_genomic_to_iso(&gff, gene, &position, &isopos);
  
  PROTECT(result=R_splicing_matrix_int_to_SEXP(&isopos));
  
  splicing_matrix_int_destroy(&isopos);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_matchIso_paired(SEXP pgff, SEXP pgene, SEXP pposition,
				SEXP pcigar, SEXP preadlength, 
				SEXP pfragmentprob, SEXP pfragmentstart, 
				SEXP pnormalMean, SEXP pnormalVar,
				SEXP pnumDevs) {
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t position;
  size_t i, noreads=GET_LENGTH(pposition);
  const char **cigarstr;
  int readlength=INTEGER(preadlength)[0];
  splicing_vector_t fragmentprob;
  int fragmentstart=INTEGER(pfragmentstart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  splicing_matrix_t result;
  splicing_matrix_int_t fragmentLength;
  SEXP rresult;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  cigarstr = (const char**) R_alloc(noreads, sizeof(char*));
  for (i=0; i<noreads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }
  if (!isNull(pfragmentprob)) {
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  splicing_matrix_init(&result, 0, 0);
  splicing_matrix_int_init(&fragmentLength, 0, 0);
  
  splicing_matchIso_paired(&gff, gene, &position, cigarstr, readlength, 
			   isNull(pfragmentprob) ? 0 : &fragmentprob, 
			   fragmentstart, normalMean, normalVar, numDevs,
			   &result, &fragmentLength);

  PROTECT(rresult=NEW_LIST(2));
  SET_VECTOR_ELT(rresult, 0, R_splicing_matrix_to_SEXP(&result));
  SET_VECTOR_ELT(rresult, 1, R_splicing_matrix_int_to_SEXP(&fragmentLength));

  splicing_matrix_destroy(&result);
  splicing_matrix_int_destroy(&fragmentLength);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return rresult;
}

SEXP R_splicing_paired_assignment_matrix(SEXP pgff, SEXP pgene, 
					 SEXP preadlength, 
					 SEXP pfragmentprob,
					 SEXP pfragmentstart, 
					 SEXP pnormalMean, SEXP pnormalVar,
					 SEXP pnumDevs) {
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readlength=INTEGER(preadlength)[0];
  splicing_vector_t fragmentprob;
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];  
  int fragmentstart=INTEGER(pfragmentstart)[0];
  splicing_matrix_t matrix;
  SEXP result;

  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  if (!isNull(pfragmentprob)) {
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  splicing_matrix_init(&matrix, 0, 0);
  
  splicing_paired_assignment_matrix(&gff, gene, readlength, 
				    isNull(pfragmentprob) ? 0 : &fragmentprob,
				    fragmentstart, normalMean, normalVar,
				    numDevs, &matrix);
  
  PROTECT(result=R_splicing_matrix_to_SEXP(&matrix));
  
  splicing_matrix_destroy(&matrix);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_splicing_normal_fragment(SEXP pnormalMean, SEXP pnormalVar, 
				SEXP pnumDevs, SEXP pminLength) { 
  SEXP result, names;
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  int minLength=INTEGER(pminLength)[0];
  splicing_vector_t fragmentProb;
  int fragmentStart;
  
  R_splicing_begin();
  
  splicing_vector_init(&fragmentProb, 0);
  splicing_normal_fragment(normalMean, normalVar, numDevs, minLength, 
			   &fragmentProb, &fragmentStart);

  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_to_SEXP(&fragmentProb));
  splicing_vector_destroy(&fragmentProb);
  SET_VECTOR_ELT(result, 1, ScalarInteger(fragmentStart));
  PROTECT(names=NEW_CHARACTER(2));
  SET_STRING_ELT(names, 0, mkChar("fragmentProb"));
  SET_STRING_ELT(names, 1, mkChar("fragmentStart"));
  SET_NAMES(result, names);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_gene_complexity(SEXP pgff, SEXP pgene, 
				SEXP preadlength, SEXP ptype, 
				SEXP pnorm, SEXP ppaired, SEXP pfragmentprob,
				SEXP pfragmentstart, SEXP pnormalmean,
				SEXP pnormalvar, SEXP pnumdevs) {
  SEXP result;
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readlength=INTEGER(preadlength)[0];
  splicing_complexity_t type=INTEGER(ptype)[0];
  splicing_norm_t norm=INTEGER(pnorm)[0];
  int paired=LOGICAL(ppaired)[0];
  splicing_vector_t fragmentprob;
  int fragmentstart=INTEGER(pfragmentstart)[0];
  double normalmean=REAL(pnormalmean)[0];
  double normalvar=REAL(pnormalvar)[0];
  double numdevs=REAL(pnumdevs)[0];
  double complexity;

  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgff, &gff);
  if (!isNull(pfragmentprob)) {
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  
  splicing_gene_complexity(&gff, gene, readlength, type, norm, paired,
			   isNull(pfragmentprob) ? 0 : &fragmentprob,
			   fragmentstart, normalmean, normalvar, numdevs,
			   &complexity);
  
  PROTECT(result=ScalarReal(complexity));

  R_splicing_end();

  UNPROTECT(1);  
  return result;
}

SEXP R_splicing_read_sambam(SEXP pfilename) {
  const char *filename=CHAR(STRING_ELT(pfilename, 0));
  splicing_reads_t reads;
  SEXP result;
  
  R_splicing_begin();

  splicing_reads_init(&reads);
  splicing_read_sambam(filename, SPLICING_SAMBAM_AUTO, &reads);
  PROTECT(result=R_splicing_reads_to_SEXP(&reads));
  splicing_reads_destroy(&reads);
  
  R_splicing_end();

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_read_sambam_region(SEXP pfilename, SEXP pregion) {
  const char *filename=CHAR(STRING_ELT(pfilename, 0));
  const char *region=CHAR(STRING_ELT(pregion, 0));
  splicing_reads_t reads;
  SEXP result;
  
  R_splicing_begin();
  
  splicing_reads_init(&reads);
  splicing_read_sambam_region(filename, /*indexfile=*/ 0, 
			      SPLICING_SAMBAM_AUTO, region, &reads);
  PROTECT(result=R_splicing_reads_to_SEXP(&reads));
  splicing_reads_destroy(&reads);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}
