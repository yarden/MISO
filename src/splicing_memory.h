/* -*- mode: C -*-  */

#ifndef SPLICING_MEMORY_H
#define SPLICING_MEMORY_H

#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#define splicing_Calloc(n,t)    (t*) calloc( (size_t)(n), sizeof(t) )
#define splicing_Realloc(p,n,t) (t*) realloc((void*)(p), (size_t)((n)*sizeof(t)))
#define splicing_Free(p)        (free( (void *)(p) ), (p) = NULL)

/* #ifndef SPLICING_NO_CALLOC */
/* #  define Calloc splicing_Calloc */
/* #  define Realloc splicing_Realloc */
/* #  define Free splicing_Free */
/* #endif */

int splicing_free(void *p);

__END_DECLS

#endif
