/* -*- mode: C -*-  */

#ifndef SPLICING_VECTOR_H
#define SPLICING_VECTOR_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#ifdef HAVE_STDINT_H
#  include <stdint.h>
#else
#  if HAVE_SYS_INT_TYPES_H
#    include <sys/int_types.h>    /* for Solaris */
#  endif
#endif

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Flexible vector                                    */
/* -------------------------------------------------- */

#define BASE_DOUBLE
#include "splicing_pmt.h"
#include "splicing_vector_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_DOUBLE

#define BASE_LONG
#include "splicing_pmt.h"
#include "splicing_vector_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "splicing_pmt.h"
#include "splicing_vector_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_CHAR

#define BASE_INT
#include "splicing_pmt.h"
#include "splicing_vector_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_INT

int splicing_vector_floor(const splicing_vector_t *from, splicing_vector_long_t *to);
int splicing_vector_round(const splicing_vector_t *from, splicing_vector_long_t *to);

/* These are for internal use only */
int splicing_vector_order(const splicing_vector_t* v, const splicing_vector_t *v2,
			splicing_vector_t* res, double maxval);
int splicing_vector_order1(const splicing_vector_t* v, 
			 splicing_vector_t* res, double maxval);
int splicing_vector_order2(splicing_vector_t *v);
int splicing_vector_rank(const splicing_vector_t *v, splicing_vector_t *res, 
		       long int nodes);

__END_DECLS

#endif
