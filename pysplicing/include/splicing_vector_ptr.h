/* -*- mode: C -*-  */

#ifndef SPLICING_VECTOR_PTR_H
#define SPLICING_VECTOR_PTR_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "splicing_vector.h"
#include "splicing_error.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Flexible vector, storing pointers                  */
/* -------------------------------------------------- */

/**
 * Vector, storing pointers efficiently
 * \ingroup internal
 *
 */
typedef struct s_vector_ptr {
  void** stor_begin;
  void** stor_end;
  void** end;
  splicing_finally_func_t* item_destructor;
} splicing_vector_ptr_t;

#define SPLICING_VECTOR_PTR_NULL { 0,0,0,0 }
#define SPLICING_VECTOR_PTR_INIT_FINALLY(v, size) \
  do { SPLICING_CHECK(splicing_vector_ptr_init(v, size)); \
  SPLICING_FINALLY(splicing_vector_ptr_destroy, v); } while (0)

int splicing_vector_ptr_init      (splicing_vector_ptr_t* v, long int size);
int splicing_vector_ptr_init_copy (splicing_vector_ptr_t* v, void** data, long int length);
const splicing_vector_ptr_t *splicing_vector_ptr_view (const splicing_vector_ptr_t *v,
				     void *const *data, long int length);
void splicing_vector_ptr_destroy   (splicing_vector_ptr_t* v);
void splicing_vector_ptr_free_all   (splicing_vector_ptr_t* v);
void splicing_vector_ptr_destroy_all   (splicing_vector_ptr_t* v);
int splicing_vector_ptr_reserve   (splicing_vector_ptr_t* v, long int size);
int splicing_vector_ptr_empty     (const splicing_vector_ptr_t* v);
long int splicing_vector_ptr_size      (const splicing_vector_ptr_t* v);
void splicing_vector_ptr_clear     (splicing_vector_ptr_t* v);
void splicing_vector_ptr_null      (splicing_vector_ptr_t* v);
int splicing_vector_ptr_push_back (splicing_vector_ptr_t* v, void* e);
int splicing_vector_ptr_append    (splicing_vector_ptr_t *to,
				 const splicing_vector_ptr_t *from);
void *splicing_vector_ptr_pop_back (splicing_vector_ptr_t *v);
int splicing_vector_ptr_insert(splicing_vector_ptr_t *v, long int pos, void* e);
void* splicing_vector_ptr_e         (const splicing_vector_ptr_t* v, long int pos);
void splicing_vector_ptr_set       (splicing_vector_ptr_t* v, long int pos, void* value);
int splicing_vector_ptr_resize(splicing_vector_ptr_t* v, long int newsize);
void splicing_vector_ptr_copy_to(const splicing_vector_ptr_t *v, void** to);
int splicing_vector_ptr_copy(splicing_vector_ptr_t *to, const splicing_vector_ptr_t *from);
void splicing_vector_ptr_remove(splicing_vector_ptr_t *v, long int pos);
void splicing_vector_ptr_sort(splicing_vector_ptr_t *v, int(*compar)(const void*, const void*));
int splicing_vector_ptr_index_int(splicing_vector_ptr_t *v,
				const splicing_vector_int_t *idx);

splicing_finally_func_t* splicing_vector_ptr_get_item_destructor(const splicing_vector_ptr_t *v);
splicing_finally_func_t* splicing_vector_ptr_set_item_destructor(splicing_vector_ptr_t *v,
        splicing_finally_func_t *func);

/**
 * \define SPLICING_VECTOR_PTR_SET_ITEM_DESTRUCTOR
 * \brief Sets the item destructor for this pointer vector (macro version).
 *
 * This macro is expanded to \ref splicing_vector_ptr_set_item_destructor(), the
 * only difference is that the second argument is automatically cast to an
 * \c splicing_finally_func_t*. The cast is necessary in most cases as the
 * destructor functions we use (such as \ref splicing_vector_destroy()) take a
 * pointer to some concrete splicing data type, while \c splicing_finally_func_t
 * expects \c void*
 */
#define SPLICING_VECTOR_PTR_SET_ITEM_DESTRUCTOR(v, func) \
        splicing_vector_ptr_set_item_destructor((v), (splicing_finally_func_t*)(func))

__END_DECLS

#endif
