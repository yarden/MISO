
#ifndef PYSPLICING_ERROR_H
#define PYSPLICING_ERROR_H

#include <Python.h>
#include "splicing_error.h"

PyObject* splicingmodule_handle_splicing_error(void);
void splicingmodule_splicing_warning_hook(const char *reason, 
					  const char *file,
					  int line, int splicing_errno);
void splicingmodule_splicing_error_hook(const char *reason, const char *file,
					int line, int splicing_errno);

extern PyObject* splicingmodule_InternalError;

#define SPLICING_PYCHECK(a) do {			 \
    int splicing_i_pyret=(a);				 \
    if (SPLICING_UNLIKELY(splicing_i_pyret != 0)) {	 \
      splicingmodule_handle_splicing_error();		 \
      SPLICING_FINALLY_FREE();				 \
      return 0;						 \
    } } while (0)

#endif
