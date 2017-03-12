
#ifndef PYSPLICING_H
#define PYSPLICING_H

#include <Python.h>

#include "splicing.h"
#include "splicing_memory.h"

int pysplicing_to_vector_int(PyObject *pv, splicing_vector_int_t *v);
int pysplicing_to_vector(PyObject *pv, splicing_vector_t *v);
int pysplicing_to_strvector(PyObject *pv, splicing_strvector_t *v);
int pysplicing_to_exons(PyObject *pex, splicing_vector_int_t *ex);
int pysplicing_to_isoforms(PyObject *piso, splicing_vector_int_t *iso);
int pysplicing_to_reads(PyObject *pyreads, splicing_reads_t *reads);
int pysplicing_to_replicate_reads(PyObject *pyreads,
				  splicing_replicate_reads_t *reads);

PyObject *pysplicing_from_vector(const splicing_vector_t *v);
PyObject *pysplicing_from_vector_int(const splicing_vector_int_t *v);
PyObject *pysplicing_from_matrix(const splicing_matrix_t *m);
PyObject *pysplicing_from_strvector(const splicing_strvector_t *v);
PyObject *pysplicing_from_vector_int_index(const splicing_vector_int_t *v,
					   const splicing_vector_int_t *idx);
PyObject *pysplicing_from_miso_rundata(const splicing_miso_rundata_t *data);
PyObject *pysplicing_from_vectorlist(const splicing_vector_ptr_t *v);
PyObject *pysplicing_from_vectorlist_int(const splicing_vector_ptr_t *v);
PyObject *pysplicing_from_matrixlist(const splicing_vector_ptr_t *v);

void pysplicing_init_rng(PyObject* splicing_module);

#endif
