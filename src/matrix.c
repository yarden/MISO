/* -*- mode: C -*-  */

#include "splicing_matrix.h"

#define BASE_DOUBLE
#include "splicing_pmt.h"
#include "matrix.pmt"
#include "splicing_pmt_off.h"
#undef BASE_DOUBLE

#define BASE_LONG
#include "splicing_pmt.h"
#include "matrix.pmt"
#include "splicing_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "splicing_pmt.h"
#include "matrix.pmt"
#include "splicing_pmt_off.h"
#undef BASE_CHAR

#define BASE_INT
#include "splicing_pmt.h"
#include "matrix.pmt"
#include "splicing_pmt_off.h"
#undef BASE_INT

