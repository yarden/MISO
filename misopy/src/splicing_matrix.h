/* -*- mode: C -*-  */

#ifndef SPLICING_MATRIX_H
#define SPLICING_MATRIX_H

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

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Matrix, very similar to vector                     */
/* -------------------------------------------------- */

#define BASE_DOUBLE
#include "splicing_pmt.h"
#include "splicing_matrix_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_DOUBLE

#define BASE_LONG
#include "splicing_pmt.h"
#include "splicing_matrix_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "splicing_pmt.h"
#include "splicing_matrix_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_CHAR

#define BASE_INT
#include "splicing_pmt.h"
#include "splicing_matrix_pmt.h"
#include "splicing_pmt_off.h"
#undef BASE_INT

#define SPLICING_MATRIX_NULL { SPLICING_VECTOR_NULL, 0, 0 }
#define SPLICING_MATRIX_INIT_FINALLY(m, nr, nc) \
  do { SPLICING_CHECK(splicing_matrix_init(m, nr, nc)); \
  SPLICING_FINALLY(splicing_matrix_destroy, m); } while (0)

/**
 * \ingroup matrix
 * \define MATRIX
 * \brief Accessing an element of a matrix.
 *
 * Note that there are no range checks right now. 
 * This functionality might be redefines as a proper function later. 
 * \param m The matrix object.
 * \param i The index of the row, starting with zero.
 * \param j The index of the column, starting with zero.
 *
 * Time complexity: O(1).
 */
#define MATRIX(m,i,j) ((m).data.stor_begin[(m).nrow*(j)+(i)])

__END_DECLS

#endif
