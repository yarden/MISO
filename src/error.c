/* -*- mode: C -*-  */

#include "splicing_error.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static splicing_error_handler_t *splicing_i_error_handler=0;

static char *splicing_i_error_strings[]=
  { /*  0 */ "No error",
    /*  1 */ "Failed",
    /*  2 */ "Out of memory",
    /*  3 */ "Parse error",
    /*  4 */ "Invalid value",
    /*  5 */ "Already exists",
    /*  6 */ "Invalid edge vector",
    /*  7 */ "Invalid vertex id",
    /*  8 */ "Non-square matrix",
    /*  9 */ "Invalid mode",
    /* 10 */ "File operation error",
    /* 11 */ "Unfold infinite iterator",
    /* 12 */ "Unimplemented function call",
    /* 13 */ "Interrupted",
    /* 14 */ "Numeric procedure did not converge",
    /* 15 */ "Matrix-vector product failed",
    /* 16 */ "N must be positive", 
    /* 17 */ "NEV must be positive",
    /* 18 */ "NCV must be bigger",
    /* 19 */ "Maximum number of iterations should be positive",
    /* 20 */ "Invalid WHICH parameter",
    /* 21 */ "Invalid BMAT parameter",
    /* 22 */ "WORKL is too small",
    /* 23 */ "LAPACK error in tridiagonal eigenvalue calculation",
    /* 24 */ "Starting vector is zero",
    /* 25 */ "MODE is invalid",
    /* 26 */ "MODE and BMAT are not compatible",
    /* 27 */ "ISHIFT must be 0 or 1",
    /* 28 */ "NEV and WHICH='BE' are incompatible",
    /* 29 */ "Could not build an Arnoldi factorization",
    /* 30 */ "No eigenvalues to sufficient accuracy",
    /* 31 */ "HOWMNY is invalid",
    /* 32 */ "HOWMNY='S' is not implemented",
    /* 33 */ "Different number of converged Ritz values",
    /* 34 */ "Error from calculation of a real Schur form",
    /* 35 */ "LAPACK (dtrevc) error for calculating eigenvectors",
    /* 36 */ "Unknown ARPACK error",
    /* 37 */ "Negative loop detected while calculating shortest paths",
    /* 38 */ "Internal error, likely a bug in splicing",
    /* 39 */ "Maximum number of iterations reached",
    /* 40 */ "No shifts could be applied during a cycle of the "
             "Implicitly restarted Arnoldi iteration. One possibility "
             "is to increase the size of NCV relative to NEV",
    /* 41 */ "The Schur form computed by LAPACK routine dlahqr "
             "could not be reordered by LAPACK routine dtrsen.",
    /* 42 */ "Big integer division by zero",
    /* 43 */ "GLPK Error, GLP_EBOUND",
    /* 44 */ "GLPK Error, GLP_EROOT",
    /* 45 */ "GLPK Error, GLP_ENOPFS",
    /* 46 */ "GLPK Error, GLP_ENODFS",
    /* 47 */ "GLPK Error, GLP_EFAIL",
    /* 48 */ "GLPK Error, GLP_EMIPGAP",
    /* 49 */ "GLPK Error, GLP_ETMLIM",
    /* 50 */ "GLPK Error, GLP_STOP",
    /* 51 */ "Internal attribute handler error",
    /* 52 */ "Unimplemented attribute combination for this type",
    /* 53 */ "LAPACK call resulted an error"
};

const char* splicing_strerror(const int splicing_errno) {
  return splicing_i_error_strings[splicing_errno];
}

int splicing_error(const char *reason, const char *file, int line,
		 int splicing_errno) {

  if (splicing_i_error_handler) {
    splicing_i_error_handler(reason, file, line, splicing_errno);
  }  else {
    splicing_error_handler_abort(reason, file, line, splicing_errno);
  }
  return splicing_errno;
}

void splicing_error_handler_abort (const char *reason, const char *file,
				 int line, int splicing_errno) {
  fprintf(stderr, "Error at %s:%i :%s, %s\n", file, line, reason,
	  splicing_strerror(splicing_errno));
  abort();
}

void splicing_error_handler_ignore (const char *reason, const char *file,
				  int line, int splicing_errno) {
  SPLICING_FINALLY_FREE();
}

void splicing_error_handler_printignore (const char *reason, const char *file,
				       int line, int splicing_errno) {
  SPLICING_FINALLY_FREE();
  fprintf(stderr, "Error at %s:%i :%s, %s\n", file, line, reason,
	  splicing_strerror(splicing_errno));
}

splicing_error_handler_t *
splicing_set_error_handler (splicing_error_handler_t * new_handler)
{
  splicing_error_handler_t * previous_handler = splicing_i_error_handler;
  splicing_i_error_handler = new_handler;
  return previous_handler;
}

struct splicing_i_protectedPtr splicing_i_finally_stack[100];

/*
 * Adds another element to the free list
 */

void SPLICING_FINALLY_REAL(void (*func)(void*), void* ptr) {
  int no=splicing_i_finally_stack[0].all;
  assert (no<100);
  assert (no>=0);
  splicing_i_finally_stack[no].ptr=ptr;
  splicing_i_finally_stack[no].func=func;
  splicing_i_finally_stack[0].all ++;
  /* printf("--> Finally stack contains now %d elements\n", splicing_i_finally_stack[0].all); */
}

void SPLICING_FINALLY_CLEAN(int minus) { 
  splicing_i_finally_stack[0].all -= minus;
  if (splicing_i_finally_stack[0].all < 0) {
    fprintf(stderr, "corrupt finally stack, popping %d elements when only %d left\n", minus, splicing_i_finally_stack[0].all+minus);
    splicing_i_finally_stack[0].all = 0;
  }
  /* printf("<-- Finally stack contains now %d elements\n", splicing_i_finally_stack[0].all); */
}

void SPLICING_FINALLY_FREE(void) {
  int p;
/*   printf("[X] Finally stack will be cleaned (contained %d elements)\n", splicing_i_finally_stack[0].all);  */
  for (p=splicing_i_finally_stack[0].all-1; p>=0; p--) {
    splicing_i_finally_stack[p].func(splicing_i_finally_stack[p].ptr);
  }
  splicing_i_finally_stack[0].all=0;
}

int SPLICING_FINALLY_STACK_SIZE(void) {
  return splicing_i_finally_stack[0].all;
}

static splicing_warning_handler_t *splicing_i_warning_handler=0;

void splicing_warning_handler_ignore (const char *reason, const char *file,
				   int line, int splicing_errno) {
}

void splicing_warning_handler_print (const char *reason, const char *file,
				   int line, int splicing_errno) {
  fprintf(stderr, "Warning: %s in file %s, line %i\n", reason, file, line);
}

int splicing_warning(const char *reason, const char *file, int line,
		   int splicing_errno) {

  if (splicing_i_warning_handler) {
    splicing_i_warning_handler(reason, file, line, splicing_errno);
  }  else {
    splicing_warning_handler_print(reason, file, line, splicing_errno);
  }
  return splicing_errno;
}

splicing_warning_handler_t *
splicing_set_warning_handler (splicing_warning_handler_t * new_handler)
{
  splicing_warning_handler_t * previous_handler = splicing_i_warning_handler;
  splicing_i_warning_handler = new_handler;
  return previous_handler;
}
