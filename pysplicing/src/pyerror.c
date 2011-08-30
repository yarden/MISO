
#include "pyerror.h"
#include "pysplicing.h"

PyObject* splicingmodule_InternalError;

PyObject* splicingmodule_handle_splicing_error() 
{
  if (!PyErr_Occurred()) {
    PyErr_SetString(splicingmodule_InternalError,
		    "Internal splicing error. Please contact the author!");
  }

  return NULL;
}

void splicingmodule_splicing_warning_hook(const char *reason, 
					  const char *file,
					  int line, int splicing_errno) {
  char buf[4096];
  snprintf(buf, sizeof(buf)/sizeof(char)-1, 
	   "%s at %s:%i", reason, file, line);
  PyErr_Warn(PyExc_RuntimeWarning, buf);
}

void splicingmodule_splicing_error_hook(const char *reason, const char *file,
					int line, int splicing_errno) {
  char buf[4096];
  PyObject *exc = splicingmodule_InternalError;

  if (splicing_errno == SPLICING_UNIMPLEMENTED)
      exc = PyExc_NotImplementedError;

  if (splicing_errno == SPLICING_ENOMEM)
      exc = PyExc_MemoryError;

  snprintf(buf, sizeof(buf)/sizeof(char)-1, 
	   "Error at %s:%i: %s, %s", file, line, reason,
	   splicing_strerror(splicing_errno));
  SPLICING_FINALLY_FREE();

  /* make sure we are not masking already thrown exceptions */
  if (!PyErr_Occurred())
    PyErr_SetString(exc, buf);
}
