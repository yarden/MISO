
#include "pysplicing.h"

int pysplicing_to_vector_int(PyObject *pv, splicing_vector_int_t *v) {
  int i, n;
  
  if (!PyTuple_Check(pv)) { 
    PyErr_SetString(PyExc_TypeError, "Need a tuple");
    return 1;
  }
  
  n=PyTuple_Size(pv);
  splicing_vector_int_init(v, n);
  for (i=0; i<n; i++) {
    PyObject *it=PyTuple_GetItem(pv, i);
    VECTOR(*v)[i]=PyInt_AsLong(it);
  }
  
  return 0;
}

int pysplicing_to_vector(PyObject *pv, splicing_vector_t *v) {
  int i, n;
  
  if (!PyTuple_Check(pv)) {
    PyErr_SetString(PyExc_TypeError, "Need a tuple");
    return 1;
  }
  
  n=PyTuple_Size(pv);
  splicing_vector_init(v, n);
  for (i=0; i<n; i++) {
    PyObject *it=PyTuple_GetItem(pv, i);
    VECTOR(*v)[i]=PyFloat_AsDouble(it);
  }
  
  return 0;
}

int pysplicing_to_strvector(PyObject *pv, splicing_strvector_t *v) {
  int i, n;

  if (!PyTuple_Check(pv)) {
    PyErr_SetString(PyExc_TypeError, "Need a tuple");
    return 1;
  }
  
  n=PyTuple_Size(pv);
  splicing_strvector_init(v, 0);
  splicing_strvector_reserve(v, n);
  for (i=0; i<n; i++) {
    PyObject *it=PyTuple_GetItem(pv, i);
    splicing_strvector_append(v,PyString_AsString(it));
  }
  
  return 0;
}

int pysplicing_to_exons(PyObject *pex, splicing_vector_int_t *ex) {
  int i, p, noexons=PyTuple_Size(pex);
  splicing_vector_int_init(ex, noexons*2);
  for (i=0, p=0; i<noexons; i++) {
    PyObject *it=PyTuple_GetItem(pex, i);
    VECTOR(*ex)[p++] = PyInt_AsLong(PyTuple_GetItem(it, 0));
    VECTOR(*ex)[p++] = PyInt_AsLong(PyTuple_GetItem(it, 1));
  }
  return 0;
}

int pysplicing_to_isoforms(PyObject *piso, splicing_vector_int_t *iso) {
  int i, p, veclen, noiso=PyTuple_Size(piso);
  
  for (i=0, veclen=0; i<noiso; i++) {
    veclen += PyTuple_Size(PyTuple_GetItem(piso, i));
  }
  
  splicing_vector_int_init(iso, veclen+noiso);
  
  for (i=0, p=0; i<noiso; i++) {
    PyObject *it=PyTuple_GetItem(piso, i);
    int j, ilen=PyTuple_Size(it);
    for (j=0; j<ilen; j++) {
      VECTOR(*iso)[p++] = PyInt_AsLong(PyTuple_GetItem(it, j));
    }
    VECTOR(*iso)[p++] = -1;
  }
  
  return 0;
}

PyObject *pysplicing_from_vector(const splicing_vector_t *v) {
  int i, n=splicing_vector_size(v);
  PyObject *o=PyTuple_New(n);
  for (i=0; i<n; i++) {
    PyObject *it=PyFloat_FromDouble(VECTOR(*v)[i]);
    PyTuple_SetItem(o, i, it);
  }
  return o;
}

PyObject *pysplicing_from_matrix(const splicing_matrix_t *m) {
  int i, j;
  int nrow=splicing_matrix_nrow(m);
  int ncol=splicing_matrix_ncol(m);
  PyObject *o=PyTuple_New(nrow);
  for (i=0; i<nrow; i++) {
    PyObject *r=PyTuple_New(ncol);
    for (j=0; j<ncol; j++) {
      PyObject *it=PyFloat_FromDouble(MATRIX(*m, i, j));
      PyTuple_SetItem(r, j, it);
    }
    PyTuple_SetItem(o, i, r);
  }
  
  return o;
}

PyObject *pysplicing_from_matrix_transposed(const splicing_matrix_t *m) {
  int i, j;
  int nrow=splicing_matrix_nrow(m);
  int ncol=splicing_matrix_ncol(m);
  PyObject *o=PyTuple_New(ncol);
  for (i=0; i<ncol; i++) {
    PyObject *r=PyTuple_New(nrow);
    for (j=0; j<nrow; j++) {
      PyObject *it=PyFloat_FromDouble(MATRIX(*m, j, i));
      PyTuple_SetItem(r, j, it);
    }
    PyTuple_SetItem(o, i, r);
  }
  
  return o;
}

PyObject *pysplicing_from_vector_int(const splicing_vector_int_t *v) {
  int i, n=splicing_vector_int_size(v);
  PyObject *o=PyTuple_New(n);
  for (i=0; i<n; i++) {
    PyObject *it=PyInt_FromLong(VECTOR(*v)[i]);
    PyTuple_SetItem(o, i, it);
  }
  return o;
}


PyObject *pysplicing_from_strvector(const splicing_strvector_t *v) {
  int i, n=splicing_strvector_size(v);
  PyObject *o=PyTuple_New(n);
  for (i=0; i<n; i++) {
    PyObject *it=PyString_FromString(splicing_strvector_get(v, i));
    PyTuple_SetItem(o, i, it);
  }
  return o;
}

PyObject *pysplicing_from_vector_int_index(const splicing_vector_int_t *v,
					   const splicing_vector_int_t *idx) {
  int i, n=splicing_vector_int_size(idx);
  int vn=splicing_vector_int_size(v);
  PyObject *o=PyList_New(n);
  for (i=0; i<n; i++) {
    int j=VECTOR(*idx)[i];
    int x, k= i<n-1 ? VECTOR(*idx)[i+1] : vn;
    PyObject *lo=PyList_New(k-j);
    for (x=0; j<k; j++, x++) {
      PyObject *it=PyInt_FromLong(VECTOR(*v)[j]);
      PyList_SetItem(lo, x, it);
    }
    PyList_SetItem(o, i, lo);
  }
  return o;
}

PyObject *pysplicing_from_miso_rundata(const splicing_miso_rundata_t *data) {
  PyObject *o=PyTuple_New(6);
  PyTuple_SetItem(o, 0, PyInt_FromLong(data->noIso));
  PyTuple_SetItem(o, 1, PyInt_FromLong(data->noIters));
  PyTuple_SetItem(o, 2, PyInt_FromLong(data->noBurnIn));
  PyTuple_SetItem(o, 3, PyInt_FromLong(data->noLag));
  PyTuple_SetItem(o, 4, PyInt_FromLong(data->noAccepted));
  PyTuple_SetItem(o, 5, PyInt_FromLong(data->noRejected));
  return o;
}
