/* -*- mode: C -*-  */

#include "splicing_vector_ptr.h"
#include "splicing_memory.h"
#include "splicing_random.h"
#include "splicing_error.h"

#include <assert.h>
#include <string.h> 		/* memcpy & co. */
#include <stdlib.h>

int splicing_vector_ptr_init      (splicing_vector_ptr_t* v, int long size) {
        long int alloc_size= size > 0 ? size : 1;
	assert(v != NULL);
	if (size < 0) { size=0; }
	v->stor_begin=splicing_Calloc(alloc_size, void*);
	if (v->stor_begin==0) {
	  SPLICING_ERROR("vector ptr init failed", SPLICING_ENOMEM);
	}
	v->stor_end=v->stor_begin + alloc_size;
	v->end=v->stor_begin+size;
	v->item_destructor=0;

	return 0;
}

const splicing_vector_ptr_t *splicing_vector_ptr_view (const splicing_vector_ptr_t *v, void *const *data,
				     long int length) {
  splicing_vector_ptr_t *v2=(splicing_vector_ptr_t*) v;
  v2->stor_begin=(void **)data;
  v2->stor_end=(void**)data+length;
  v2->end=v2->stor_end;
  v2->item_destructor=0;
  return v;
}


void splicing_vector_ptr_destroy   (splicing_vector_ptr_t* v) {
  assert(v != 0);
  if (v->stor_begin != 0) {
    splicing_Free(v->stor_begin);
    v->stor_begin = NULL;
  }
}

void splicing_i_vector_ptr_call_item_destructor_all(splicing_vector_ptr_t* v) {
  void **ptr;

  if (v->item_destructor != 0) {
    for (ptr=v->stor_begin; ptr<v->end; ptr++) {
      if (*ptr != 0)
        v->item_destructor(*ptr);
    }
  }
}

void splicing_vector_ptr_free_all   (splicing_vector_ptr_t* v) {
  void **ptr;
  assert(v != 0);
  assert(v->stor_begin != 0);

  splicing_i_vector_ptr_call_item_destructor_all(v);
  for (ptr=v->stor_begin; ptr<v->end; ptr++) {
    splicing_Free(*ptr);
  }
}


void splicing_vector_ptr_destroy_all   (splicing_vector_ptr_t* v) {
  assert(v != 0);
  assert(v->stor_begin != 0);
  splicing_vector_ptr_free_all(v);
  splicing_vector_ptr_set_item_destructor(v, 0);
  splicing_vector_ptr_destroy(v);
}


int splicing_vector_ptr_reserve   (splicing_vector_ptr_t* v, long int size) {
	long int actual_size=splicing_vector_ptr_size(v);
	void **tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);

	if (size <= splicing_vector_ptr_size(v)) { return 0; }

	tmp=splicing_Realloc(v->stor_begin, (size_t) size, void*);
	if (tmp==0) {
	  SPLICING_ERROR("vector ptr reserve failed", SPLICING_ENOMEM);
	}
	v->stor_begin=tmp;
	v->stor_end=v->stor_begin + size;
	v->end=v->stor_begin+actual_size;

	return 0;
}


int splicing_vector_ptr_empty     (const splicing_vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return v->stor_begin == v->end;
}


long int splicing_vector_ptr_size      (const splicing_vector_ptr_t* v) {
	assert(v != NULL);
/* 	assert(v->stor_begin != NULL);		 */ /* TODO */
	return v->end - v->stor_begin;
}


void splicing_vector_ptr_clear     (splicing_vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	splicing_i_vector_ptr_call_item_destructor_all(v);
	v->end = v->stor_begin;
}


int splicing_vector_ptr_push_back (splicing_vector_ptr_t* v, void* e) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);

	/* full, allocate more storage */
	if (v->stor_end == v->end) {
		long int new_size = splicing_vector_ptr_size(v) * 2;
		if (new_size == 0) { new_size = 1; }
		SPLICING_CHECK(splicing_vector_ptr_reserve(v, new_size));
	}

	*(v->end) = e;
	v->end += 1;

	return 0;
}

void *splicing_vector_ptr_pop_back (splicing_vector_ptr_t *v) {
	void *tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	assert(v->stor_begin != v->end);
	tmp=*(v->end);
	v->end -= 1;

	return tmp;
}

int splicing_vector_ptr_insert(splicing_vector_ptr_t* v, long int pos, void* e) {
  long int size = splicing_vector_ptr_size(v);
  SPLICING_CHECK(splicing_vector_ptr_resize(v, size+1));
  if (pos<size) {
    memmove(v->stor_begin+pos+1, v->stor_begin+pos,
	    sizeof(void*) * (size_t) (size-pos));
  }
  v->stor_begin[pos] = e;
  return 0;
}


void* splicing_vector_ptr_e         (const splicing_vector_ptr_t* v, long int pos) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return * (v->stor_begin + pos);
}


void splicing_vector_ptr_set       (splicing_vector_ptr_t* v, long int pos, void* value) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	*(v->stor_begin + pos) = value;
}


void splicing_vector_ptr_null      (splicing_vector_ptr_t* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	if (splicing_vector_ptr_size(v)>0) {
	  memset(v->stor_begin, 0, sizeof(void*) *
		 (size_t) splicing_vector_ptr_size(v));
	}
}


int splicing_vector_ptr_resize(splicing_vector_ptr_t* v, long int newsize) {
  SPLICING_CHECK(splicing_vector_ptr_reserve(v, newsize));
  v->end = v->stor_begin+newsize;
  return 0;
}


int splicing_vector_ptr_init_copy(splicing_vector_ptr_t *v, void* *data, long int length) {
  v->stor_begin=splicing_Calloc(length, void*);
  if (v->stor_begin==0) {
    SPLICING_ERROR("cannot init ptr vector from array", SPLICING_ENOMEM);
  }
  v->stor_end=v->stor_begin+length;
  v->end=v->stor_end;
  v->item_destructor=0;
  memcpy(v->stor_begin, data, (size_t) length * sizeof(void*));

  return 0;
}


void splicing_vector_ptr_copy_to(const splicing_vector_ptr_t *v, void** to) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  if (v->end != v->stor_begin) {
    memcpy(to, v->stor_begin, sizeof(void*) *
	   (size_t) (v->end - v->stor_begin));
  }
}


int splicing_vector_ptr_copy(splicing_vector_ptr_t *to, const splicing_vector_ptr_t *from) {
  assert(from != NULL);
/*   assert(from->stor_begin != NULL); */ /* TODO */
  to->stor_begin=splicing_Calloc(splicing_vector_ptr_size(from), void*);
  if (to->stor_begin==0) {
    SPLICING_ERROR("cannot copy ptr vector", SPLICING_ENOMEM);
  }
  to->stor_end=to->stor_begin+splicing_vector_ptr_size(from);
  to->end=to->stor_end;
  to->item_destructor=from->item_destructor;
  memcpy(to->stor_begin, from->stor_begin,
	 (size_t) splicing_vector_ptr_size(from)*sizeof(void*));

  return 0;
}


void splicing_vector_ptr_remove(splicing_vector_ptr_t *v, long int pos) {
  assert(v != NULL);
  assert(v->stor_begin != NULL);
  if (pos+1<splicing_vector_ptr_size(v)) { /* TOOD: why is this needed */
    memmove(v->stor_begin+pos, v->stor_begin+pos+1,
	    sizeof(void*) * (size_t) (splicing_vector_ptr_size(v)-pos-1));
  }
  v->end--;
}


void splicing_vector_ptr_sort(splicing_vector_ptr_t *v, int (*compar)(const void*, const void*)) {
  qsort(v->stor_begin, (size_t) splicing_vector_ptr_size(v), sizeof(void*),
	compar);
}

int splicing_vector_ptr_index_int(splicing_vector_ptr_t *v,
				const splicing_vector_int_t *idx) {
  void **tmp;
  int i, n=splicing_vector_int_size(idx);

  tmp=splicing_Calloc(n, void*);
  if (!tmp) { SPLICING_ERROR("Cannot index pointer vector", SPLICING_ENOMEM); }

  for (i=0; i<n; i++) { tmp[i] = VECTOR(*v)[ VECTOR(*idx)[i] ]; }

  splicing_Free(v->stor_begin);
  v->stor_begin = tmp;
  v->stor_end = v->end = tmp + n;

  return 0;
}

int splicing_vector_ptr_append    (splicing_vector_ptr_t *to,
				 const splicing_vector_ptr_t *from) {
  long int origsize=splicing_vector_ptr_size(to);
  long int othersize=splicing_vector_ptr_size(from);
  long int i;

  SPLICING_CHECK(splicing_vector_ptr_resize(to, origsize+othersize));
  for (i=0; i<othersize; i++, origsize++) {
    to->stor_begin[origsize]=from->stor_begin[i];
  }

  return 0;
}


splicing_finally_func_t* splicing_vector_ptr_set_item_destructor(
        splicing_vector_ptr_t *v, splicing_finally_func_t *func) {
  splicing_finally_func_t* result = v->item_destructor;

  v->item_destructor = func;

  return result;
}


splicing_finally_func_t* splicing_vector_ptr_get_item_destructor(const splicing_vector_ptr_t *v) {
  assert(v != 0);
  return v->item_destructor;
}
