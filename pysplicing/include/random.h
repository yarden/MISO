
#ifndef PYSPLICING_RANDOM
#define PYSPLICING_RANDOM

#include <Python.h>

void splicingmodule_init_rng(PyObject*);
PyObject* splicing_rng_Python_set_generator(PyObject* self, PyObject* object);

#endif
