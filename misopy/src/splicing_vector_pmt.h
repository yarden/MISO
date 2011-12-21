/* -*- mode: C -*-  */

/** 
 * Vector, dealing with arrays efficiently.
 * \ingroup types
 */

#include <stdio.h>

typedef struct TYPE(splicing_vector) {
  BASE* stor_begin;
  BASE* stor_end;
  BASE* end;
} TYPE(splicing_vector);

#ifndef SPLICING_VECTOR_NULL
#define SPLICING_VECTOR_NULL { 0,0,0 }
#endif

#ifndef SPLICING_VECTOR_INIT_FINALLY
#define SPLICING_VECTOR_INIT_FINALLY(v, size) \
  do { SPLICING_CHECK(splicing_vector_init(v, size)); \
  SPLICING_FINALLY(splicing_vector_destroy, v); } while (0)
#endif

#ifndef SPLICING_VECTOR_INT_INIT_FINALLY
#define SPLICING_VECTOR_INT_INIT_FINALLY(v, size) \
  do { SPLICING_CHECK(splicing_vector_int_init(v, size)); \
  SPLICING_FINALLY(splicing_vector_int_destroy, v); } while (0)
#endif

/*--------------------*/
/* Allocation         */
/*--------------------*/
 
int FUNCTION(splicing_vector,init)(TYPE(splicing_vector)* v, long int size);
int FUNCTION(splicing_vector,init_copy)(TYPE(splicing_vector)* v, 
				       BASE* data, long int length);
int FUNCTION(splicing_vector,init_seq)(TYPE(splicing_vector)*v, BASE from, BASE to);
int FUNCTION(splicing_vector,copy)(TYPE(splicing_vector) *to, 
				 const TYPE(splicing_vector) *from);
void FUNCTION(splicing_vector,destroy)(TYPE(splicing_vector)* v);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

#ifndef VECTOR
/**
 * \ingroup vector
 * \define VECTOR
 * \brief Accessing an element of a vector.
 * 
 * Usage: 
 * \verbatim VECTOR(v)[0] \endverbatim 
 * to access the first element of the vector, you can also use this in
 * assignments, like: 
 * \verbatim VECTOR(v)[10]=5; \endverbatim
 *
 * Note that there are no range checks right now.
 * This functionality might be redefined later as a real function
 * instead of a <code>#define</code>. 
 * \param v The vector object.
 * 
 * Time complexity: O(1).
 */
#define VECTOR(v) ((v).stor_begin) 
#endif

BASE FUNCTION(splicing_vector,e)(const TYPE(splicing_vector)* v, long int pos);
BASE* FUNCTION(splicing_vector,e_ptr)(const TYPE(splicing_vector)* v, long int pos);
void FUNCTION(splicing_vector,set)(TYPE(splicing_vector)* v, long int pos, BASE value);
BASE FUNCTION(splicing_vector,tail)(const TYPE(splicing_vector) *v);

/*-----------------------*/
/* Initializing elements */
/*-----------------------*/

void FUNCTION(splicing_vector,null)(TYPE(splicing_vector)* v);
void FUNCTION(splicing_vector,fill)(TYPE(splicing_vector)* v, BASE e);

/*-----------------------*/
/* Vector views          */
/*-----------------------*/

const TYPE(splicing_vector) *FUNCTION(splicing_vector,view)(const TYPE(splicing_vector) *v,
							const BASE *data, 
							long int length);

/*-----------------------*/
/* Copying vectors       */
/*-----------------------*/

void FUNCTION(splicing_vector,copy_to)(const TYPE(splicing_vector) *v, BASE* to);
int FUNCTION(splicing_vector,update)(TYPE(splicing_vector) *to, 
				   const TYPE(splicing_vector) *from);
int FUNCTION(splicing_vector,append)(TYPE(splicing_vector) *to, 
				   const TYPE(splicing_vector) *from);
int FUNCTION(splicing_vector,swap)(TYPE(splicing_vector) *v1, TYPE(splicing_vector) *v2);

/*-----------------------*/
/* Exchanging elements   */
/*-----------------------*/

int FUNCTION(splicing_vector,swap_elements)(TYPE(splicing_vector) *v,
					  long int i, long int j);
int FUNCTION(splicing_vector,reverse)(TYPE(splicing_vector) *v);
int FUNCTION(splicing_vector,shuffle)(TYPE(splicing_vector) *v);

/*-----------------------*/
/* Vector operations     */
/*-----------------------*/

void FUNCTION(splicing_vector,add_constant)(TYPE(splicing_vector) *v, BASE plus);
void FUNCTION(splicing_vector,scale)(TYPE(splicing_vector) *v, BASE by);
int FUNCTION(splicing_vector,add)(TYPE(splicing_vector) *v1, 
				const TYPE(splicing_vector) *v2);
int FUNCTION(splicing_vector,sub)(TYPE(splicing_vector) *v1, 
				const TYPE(splicing_vector) *v2);
int FUNCTION(splicing_vector,mul)(TYPE(splicing_vector) *v1, 
				const TYPE(splicing_vector) *v2);
int FUNCTION(splicing_vector,div)(TYPE(splicing_vector) *v1, 
				const TYPE(splicing_vector) *v2);

/*------------------------------*/
/* Finding minimum and maximum  */
/*------------------------------*/

BASE FUNCTION(splicing_vector,min)(const TYPE(splicing_vector)* v);
BASE FUNCTION(splicing_vector,max)(const TYPE(splicing_vector)* v);
long int FUNCTION(splicing_vector,which_min)(const TYPE(splicing_vector)* v);
long int FUNCTION(splicing_vector,which_max)(const TYPE(splicing_vector)* v);
int FUNCTION(splicing_vector,minmax)(const TYPE(splicing_vector) *v,
				   BASE *min, BASE *max);
int FUNCTION(splicing_vector,which_minmax)(const TYPE(splicing_vector) *v,
					 long int *which_min, long int *which_max);

/*-------------------*/
/* Vector properties */
/*-------------------*/

int FUNCTION(splicing_vector,empty)     (const TYPE(splicing_vector)* v);
long int FUNCTION(splicing_vector,size)      (const TYPE(splicing_vector)* v);
int FUNCTION(splicing_vector,isnull)(const TYPE(splicing_vector) *v);
BASE FUNCTION(splicing_vector,sum)(const TYPE(splicing_vector) *v);
BASE FUNCTION(splicing_vector,prod)(const TYPE(splicing_vector) *v);
int FUNCTION(splicing_vector,isininterval)(const TYPE(splicing_vector) *v, 
						   BASE low, BASE high);
int FUNCTION(splicing_vector,any_smaller)(const TYPE(splicing_vector) *v, 
						  BASE limit);
int FUNCTION(splicing_vector,is_equal)(const TYPE(splicing_vector) *lhs, 
					       const TYPE(splicing_vector) *rhs);
BASE FUNCTION(splicing_vector,maxdifference)(const TYPE(splicing_vector) *m1,
					   const TYPE(splicing_vector) *m2);

/*------------------------*/
/* Searching for elements */
/*------------------------*/

int FUNCTION(splicing_vector,contains)(const TYPE(splicing_vector) *v, BASE e);
int FUNCTION(splicing_vector,search)(const TYPE(splicing_vector) *v,
					     long int from, BASE what, 
					     long int *pos);
int FUNCTION(splicing_vector,binsearch)(const TYPE(splicing_vector) *v, 
						BASE what, long int *pos);
int FUNCTION(splicing_vector,binsearch2)(const TYPE(splicing_vector) *v,
						 BASE what);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

void FUNCTION(splicing_vector,clear)(TYPE(splicing_vector)* v);
int FUNCTION(splicing_vector,resize)(TYPE(splicing_vector)* v, long int newsize);
int FUNCTION(splicing_vector,reserve)(TYPE(splicing_vector)* v, long int size);
int FUNCTION(splicing_vector,push_back)(TYPE(splicing_vector)* v, BASE e);
BASE FUNCTION(splicing_vector,pop_back)(TYPE(splicing_vector)* v);
int FUNCTION(splicing_vector,insert)(TYPE(splicing_vector) *v, long int pos, BASE value);
void FUNCTION(splicing_vector,remove)(TYPE(splicing_vector) *v, long int elem);
void FUNCTION(splicing_vector,remove_section)(TYPE(splicing_vector) *v, 
					    long int from, long int to);

/*-----------*/
/* Sorting   */             
/*-----------*/

void FUNCTION(splicing_vector,sort)(TYPE(splicing_vector) *v);
long int FUNCTION(splicing_vector,qsort_ind)(TYPE(splicing_vector) *v, 
					   splicing_vector_t *inds, int descending);

/*-----------*/
/* Printing  */
/*-----------*/

int FUNCTION(splicing_vector,print)(const TYPE(splicing_vector) *v);
int FUNCTION(splicing_vector,fprint)(const TYPE(splicing_vector) *v, FILE *file);

#ifdef BASE_COMPLEX

int splicing_vector_complex_real(const splicing_vector_complex_t *v, 
			       splicing_vector_t *real);
int splicing_vector_complex_imag(const splicing_vector_complex_t *v, 
			       splicing_vector_t *imag);
int splicing_vector_complex_realimag(const splicing_vector_complex_t *v, 
				   splicing_vector_t *real, 
				   splicing_vector_t *imag);
int splicing_vector_complex_create(splicing_vector_complex_t *v,
				 const splicing_vector_t *real,
				 const splicing_vector_t *imag);
int splicing_vector_complex_create_polar(splicing_vector_complex_t *v,
				       const splicing_vector_t *r,
				       const splicing_vector_t *theta);

#endif

/* ----------------------------------------------------------------------------*/
/* For internal use only, may be removed, rewritten ... */
/* ----------------------------------------------------------------------------*/

int FUNCTION(splicing_vector,init_real)(TYPE(splicing_vector)*v, int no, ...);
int FUNCTION(splicing_vector,init_int)(TYPE(splicing_vector)*v, int no, ...);
int FUNCTION(splicing_vector,init_real_end)(TYPE(splicing_vector)*v, BASE endmark, ...);
int FUNCTION(splicing_vector,init_int_end)(TYPE(splicing_vector)*v, int endmark, ...);

int FUNCTION(splicing_vector,move_interval)(TYPE(splicing_vector) *v, 
					  long int begin, long int end, long int to);
int FUNCTION(splicing_vector,move_interval2)(TYPE(splicing_vector) *v, 
					  long int begin, long int end, long int to);
void FUNCTION(splicing_vector,permdelete)(TYPE(splicing_vector) *v, 
					const splicing_vector_t *index, 
					long int nremove);
void FUNCTION(splicing_vector,remove_negidx)(TYPE(splicing_vector) *v, 
					   const splicing_vector_t *neg, 
					   long int nremove);
int FUNCTION(splicing_vector,filter_smaller)(TYPE(splicing_vector) *v, BASE elem);
int FUNCTION(splicing_vector,get_interval)(const TYPE(splicing_vector) *v, 
					 TYPE(splicing_vector) *res,
					 long int from, long int to);
int FUNCTION(splicing_vector,difference_sorted)(const TYPE(splicing_vector) *v1,
  const TYPE(splicing_vector) *v2, TYPE(splicing_vector) *result);
int FUNCTION(splicing_vector,intersect_sorted)(const TYPE(splicing_vector) *v1,
  const TYPE(splicing_vector) *v2, TYPE(splicing_vector) *result);

int FUNCTION(splicing_vector,index)(const TYPE(splicing_vector) *v,
                                  TYPE(splicing_vector) *newv,
                                  const splicing_vector_t *idx);

int FUNCTION(splicing_vector,iindex)(TYPE(splicing_vector) *v,
				     const splicing_vector_t *idx);
struct splicing_vector_int_t;
int FUNCTION(splicing_vector,intiindex)(TYPE(splicing_vector) *v,
					const struct splicing_vector_int_t *idx);
