/* -*- mode: C -*-  */

#include "splicing_vector.h"
#include "splicing_memory.h"
#include "splicing.h"

#define BASE_DOUBLE
#include "splicing_pmt.h"
#include "vector.pmt"
#include "splicing_pmt_off.h"
#undef BASE_DOUBLE

#define BASE_LONG
#include "splicing_pmt.h"
#include "vector.pmt"
#include "splicing_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "splicing_pmt.h"
#include "vector.pmt"
#include "splicing_pmt_off.h"
#undef BASE_CHAR

#define BASE_INT
#include "splicing_pmt.h"
#include "vector.pmt"
#include "splicing_pmt_off.h"
#undef BASE_INT

int splicing_vector_floor(const splicing_vector_t *from, splicing_vector_long_t *to) {
  long int i, n=splicing_vector_size(from);
  
  splicing_vector_long_resize(to, n);
  for (i=0; i<n; i++) {
    VECTOR(*to)[i] = floor(VECTOR(*from)[i]);
  }
  return 0;
}

int splicing_vector_round(const splicing_vector_t *from, splicing_vector_long_t *to) {
  long int i, n=splicing_vector_size(from);
  
  splicing_vector_long_resize(to, n);
  for (i=0; i<n; i++) {
    VECTOR(*to)[i] = round(VECTOR(*from)[i]);
  }
  return 0;
}

/**
 * \ingroup vector
 * \function splicing_vector_order
 * \brief Calculate the order of the elements in a vector.
 *
 * </para><para>
 * The smallest element will have order zero, the second smallest
 * order one, etc. 
 * \param v The original \type splicing_vector_t object.
 * \param res An initialized \type splicing_vector_t object, it will be
 *    resized to match the size of \p v. The
 *    result of the computation will be stored here.
 * \param nodes Hint, the largest element in \p v.
 * \return Error code:
 *         \c SPLICING_ENOMEM: out of memory
 * 
 * Time complexity: O()
 */

int splicing_vector_order(const splicing_vector_t* v, 
			const splicing_vector_t *v2,
			splicing_vector_t* res, double nodes) {
  long int edges=splicing_vector_size(v);
  splicing_vector_t ptr;
  splicing_vector_t rad;
  long int i, j;

  assert(v!=NULL);
  assert(v->stor_begin != NULL);

  SPLICING_VECTOR_INIT_FINALLY(&ptr, nodes+1);
  SPLICING_VECTOR_INIT_FINALLY(&rad, edges);
  splicing_vector_resize(res, edges);
  
  for (i=0; i<edges; i++) {
    long int radix=v2->stor_begin[i];
    if (VECTOR(ptr)[radix]!=0) {
      VECTOR(rad)[i]=VECTOR(ptr)[radix];
    }
    VECTOR(ptr)[radix]=i+1;
  }  

  j=0;
  for (i=0; i<nodes+1; i++) {
    if (VECTOR(ptr)[i] != 0) {
      long int next=VECTOR(ptr)[i]-1;
      res->stor_begin[j++]=next;
      while (VECTOR(rad)[next] != 0) {
	next=VECTOR(rad)[next]-1;
	res->stor_begin[j++]=next;
      }
    }
  }

  splicing_vector_null(&ptr);
  splicing_vector_null(&rad);

  for (i=0; i<edges; i++) {
    long int edge=VECTOR(*res)[edges-i-1];
    long int radix=VECTOR(*v)[edge];
    if (VECTOR(ptr)[radix]!= 0) {
      VECTOR(rad)[edge]=VECTOR(ptr)[radix];
    }
    VECTOR(ptr)[radix]=edge+1;
  }
  
  j=0;
  for (i=0; i<nodes+1; i++) {
    if (VECTOR(ptr)[i] != 0) {
      long int next=VECTOR(ptr)[i]-1;
      res->stor_begin[j++]=next;
      while (VECTOR(rad)[next] != 0) {
	next=VECTOR(rad)[next]-1;
	res->stor_begin[j++]=next;
      }
    }
  } 
  
  splicing_vector_destroy(&ptr);
  splicing_vector_destroy(&rad);
  SPLICING_FINALLY_CLEAN(2);
  
  return 0;
}

int splicing_vector_order1(const splicing_vector_t* v,
			 splicing_vector_t* res, double nodes) {
  long int edges=splicing_vector_size(v);
  splicing_vector_t ptr;
  splicing_vector_t rad;
  long int i, j;

  assert(v!=NULL);
  assert(v->stor_begin != NULL);

  SPLICING_VECTOR_INIT_FINALLY(&ptr, nodes+1);
  SPLICING_VECTOR_INIT_FINALLY(&rad, edges);
  splicing_vector_resize(res, edges);
  
  for (i=0; i<edges; i++) {
    long int radix=v->stor_begin[i];
    if (VECTOR(ptr)[radix]!=0) {
      VECTOR(rad)[i]=VECTOR(ptr)[radix];
    }
    VECTOR(ptr)[radix]=i+1;
  }
  
  j=0;
  for (i=0; i<nodes+1; i++) {
    if (VECTOR(ptr)[i] != 0) {
      long int next=VECTOR(ptr)[i]-1;
      res->stor_begin[j++]=next;
      while (VECTOR(rad)[next] != 0) {
	next=VECTOR(rad)[next]-1;
	res->stor_begin[j++]=next;
      }
    }
  }
  
  splicing_vector_destroy(&ptr);
  splicing_vector_destroy(&rad);
  SPLICING_FINALLY_CLEAN(2);
  
  return 0;
}

int splicing_vector_rank(const splicing_vector_t *v, splicing_vector_t *res,
		       long int nodes) {
  
  splicing_vector_t rad;
  splicing_vector_t ptr;
  long int edges = splicing_vector_size(v);
  long int i, c=0;
  
  SPLICING_VECTOR_INIT_FINALLY(&rad, nodes);
  SPLICING_VECTOR_INIT_FINALLY(&ptr, edges);
  splicing_vector_resize(res, edges);
	       
  for (i=0; i<edges; i++) {
    long int elem=VECTOR(*v)[i];
    VECTOR(ptr)[i] = VECTOR(rad)[elem];
    VECTOR(rad)[elem] = i+1;
  }
  
  for (i=0; i<nodes; i++) {
    long int p=VECTOR(rad)[i];
    while (p != 0) {      
      VECTOR(*res)[p-1]=c++;
      p=VECTOR(ptr)[p-1];
    }
  }

  splicing_vector_destroy(&ptr);
  splicing_vector_destroy(&rad);
  SPLICING_FINALLY_CLEAN(2);
  return 0;
}
