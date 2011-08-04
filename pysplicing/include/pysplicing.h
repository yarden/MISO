
#ifndef PYSPLICING_H
#define PYSPLICING_H

#include <Python.h>

#include "splicing.h"

int pysplicing_to_vector_int(PyObject *pv, splicing_vector_int_t *v);
int pysplicing_to_vector(PyObject *pv, splicing_vector_t *v);
int pysplicing_to_strvector(PyObject *pv, splicing_strvector_t *v);
int pysplicing_to_exons(PyObject *pex, splicing_vector_int_t *ex);
int pysplicing_to_isoforms(PyObject *piso, splicing_vector_int_t *iso);

PyObject *pysplicing_from_vector(const splicing_vector_t *v);
PyObject *pysplicing_from_vector_int(const splicing_vector_int_t *v);
PyObject *pysplicing_from_matrix(const splicing_matrix_t *m);
PyObject *pysplicing_from_strvector(const splicing_strvector_t *v);
PyObject *pysplicing_from_vector_int_index(const splicing_vector_int_t *v,
					   const splicing_vector_int_t *idx);
PyObject *pysplicing_from_miso_rundata(const splicing_miso_rundata_t *data);

#endif
