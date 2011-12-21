#include "splicing.h"
#include "splicing_error.h"
#include "splicing_memory.h"

#include <string.h>
#include <ctype.h>

const char *splicing_strvector_zero="";
#define SPLICING_STRVECTOR_ZERO ((char*) splicing_strvector_zero)

int splicing_strvector_init(splicing_strvector_t *v, size_t size) {
  int i;

  v->size=size;
  v->asize=size;

  v->table = malloc(size * sizeof(char *));
  if (!v->table) { 
    SPLICING_ERROR("Cannot create string vector", SPLICING_ENOMEM);
  }
  
  for (i=0; i<size; i++) {
    v->table[i] = SPLICING_STRVECTOR_ZERO;
  }
  
  return 0;
}

int splicing_strvector_clear(splicing_strvector_t *v) {
  int i;
  for (i=0; i<v->size; i++) {
    free(v->table[i]);
  }
  v->size=0;
}

void splicing_strvector_destroy(splicing_strvector_t *v) {
  int i;
  for (i=0; i<v->size; i++) {
    free(v->table[i]);
  }
  if (v->table)  { free(v->table);  v->table=0;  }
}

size_t splicing_strvector_size(const splicing_strvector_t *v) {
  return v->size;
}

int splicing_strvector_append(splicing_strvector_t *v, const char *str) {
  if (v->asize == v->size) { 
    size_t newsize = v->size == 0 ? 1 : v->size * 2;
    SPLICING_CHECK(splicing_strvector_reserve(v, newsize));
  }

  v->table[v->size] = strdup(str);
  if (! v->table[v->size]) {
    SPLICING_ERROR("Cannot append to string vector", SPLICING_ENOMEM);
  }
  v->size++;
  
  return 0;
}

int splicing_strvector_append2(splicing_strvector_t *v, const char *str, 
			       size_t len) {
  if (v->asize == v->size) { 
    size_t newsize = v->size == 0 ? 1 : v->size * 2;
    SPLICING_CHECK(splicing_strvector_reserve(v, newsize));
  }

  v->table[v->size] = malloc( (len+1) * sizeof(char) );
  if (! v->table[v->size]) {
    SPLICING_ERROR("Cannot append to string vector", SPLICING_ENOMEM);
  }
  memcpy(v->table[v->size], str, len * sizeof(char));
  v->table[v->size][len] = '\0';
  v->size++;
  
  return 0;
}

int splicing_strvector_reserve(splicing_strvector_t *v, size_t size) {
  char **t;

  if (v->asize >= size) { return 0; }

  t=realloc(v->table, sizeof(char*) * size);

  if (!t) { 
    SPLICING_ERROR("Cannot reserve string vector", SPLICING_ENOMEM);
  }

  v->table = t;
  v->asize = size;
  
  return 0;
}

const char *splicing_strvector_get(const splicing_strvector_t *v, 
				   size_t idx) {

  if (idx >= v->size) { 
    splicing_error("String index too big", __FILE__, __LINE__, 
		   SPLICING_EINVAL);
  }
  return v->table[idx];
}

/* Return 1 if found, 0 otherwise. */

int splicing_strvector_search(const splicing_strvector_t *v, 
			      const char *key, size_t *idx) {
  int i, n=v->size;
  for (i=0; i<n; i++) {
    if (!strcmp(v->table[i], key)) { *idx = i; return 1; }
  }
  return 0;
}

int splicing_strvector_fprint(const splicing_strvector_t *v, FILE *file) {
  int i, eof=!EOF;
  for (i=0; i<v->size; i++) {
    eof = eof || fputc('"', file);
    eof = eof || fputs(v->table[i], file);
    eof = eof || fputs("\", ", file);
  }
  eof = eof || fputc('\n', file);
  if (eof == EOF) { SPLICING_ERROR("Cannot print string", SPLICING_EFILE); }
  return 0;
}

int splicing_strvector_print(const splicing_strvector_t *v) {
  return splicing_strvector_fprint(v, stdout);
}

int splicing_strvector_permute(splicing_strvector_t *v, 
			       const splicing_vector_t *idx) {

  char **copy=malloc(v->size * sizeof(char*));
  size_t i;
  
  if (!copy) { 
    SPLICING_ERROR("Cannot index string vector", SPLICING_ENOMEM); 
  }
  SPLICING_FINALLY(splicing_free, copy);
  memcpy(copy, v->table, v->size * sizeof(char*));

  for (i=0; i < v->size; i++) {
    size_t w=VECTOR(*idx)[i];
    v->table[i] = copy[w];
  }
  
  splicing_free(copy);
  SPLICING_FINALLY_CLEAN(1);

  return 0;
}

int splicing_strvector_ipermute(splicing_strvector_t *v, 
				const splicing_vector_int_t *idx) {

  char **copy=malloc(v->size * sizeof(char*));
  size_t i;
  
  if (!copy) { 
    SPLICING_ERROR("Cannot index string vector", SPLICING_ENOMEM); 
  }
  SPLICING_FINALLY(splicing_free, copy);
  memcpy(copy, v->table, v->size * sizeof(char*));

  for (i=0; i < v->size; i++) {
    size_t w=VECTOR(*idx)[i];
    v->table[i] = copy[w];
  }
  
  splicing_free(copy);
  SPLICING_FINALLY_CLEAN(1);

  return 0;
}
