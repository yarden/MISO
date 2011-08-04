
#include <Python.h>

#include <stdio.h>

#include "splicing.h"
#include "pysplicing.h"
#include "pyerror.h"
#include "random.h"
#include "splicing_memory.h"

/* TODO: exceptions */

static PyObject* pysplicing_read_gff(PyObject *self, PyObject *args) {
  const char *filename;
  FILE *input;
  splicing_gff_t *gff;

  if (!PyArg_ParseTuple(args, "s", &filename)) { return NULL; }

  input = fopen(filename, "r");
  
  gff=malloc(sizeof(splicing_gff_t));
  if (!gff) { 
    splicing_error("Cannot create GFF", __FILE__, __LINE__, SPLICING_ENOMEM); 
    splicingmodule_handle_splicing_error(); 
    return NULL; 
  } 
  SPLICING_FINALLY(splicing_free, gff);
  SPLICING_PYCHECK(splicing_gff_init(gff, 50));
  SPLICING_FINALLY(splicing_gff_destroy, gff);
  splicing_gff_read(input, gff);
  fclose(input);

  SPLICING_FINALLY_CLEAN(2);
  
  return PyCObject_FromVoidPtr(gff, splicing_gff_destroy2);
}

static PyObject* pysplicing_miso(PyObject *self, PyObject *args) {
  PyObject *gff, *readpos, *readcigar, *hyperp=0;
  int gene, readLength, noIterations=5000, noBurnIn=500, noLag=10;
  splicing_gff_t *mygff;
  splicing_strvector_t myreadcigar;
  splicing_vector_int_t myreadpos;
  splicing_vector_t myhyperp;
  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_miso_rundata_t rundata;
  PyObject *result1, *result2, *result3;
  
  if (!PyArg_ParseTuple(args, "OiOOi|iiiO", &gff, &gene, &readpos, &readcigar,
			&readLength, &noIterations, &noBurnIn, &noLag, 
			&hyperp)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);

  SPLICING_PYCHECK(splicing_matrix_init(&samples, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &samples);
  SPLICING_PYCHECK(splicing_vector_init(&logLik, 0));
  SPLICING_FINALLY(splicing_vector_destroy, &logLik);
  if (pysplicing_to_vector_int(readpos, &myreadpos)) { return NULL; }
  SPLICING_FINALLY(splicing_vector_int_destroy, &myreadpos);
  if (hyperp) { 
    if (pysplicing_to_vector(hyperp, &myhyperp)) { return NULL; }
    SPLICING_FINALLY(splicing_vector_destroy, &myhyperp);
  } else {
    size_t i, noiso;
    SPLICING_PYCHECK(splicing_gff_noiso_one(mygff, gene, &noiso));
    SPLICING_PYCHECK(splicing_vector_init(&myhyperp, noiso));
    SPLICING_FINALLY(splicing_vector_destroy, &myhyperp);
    for (i=0; i<noiso; i++) { VECTOR(myhyperp)[i] = 1.0; }
  }
  if (pysplicing_to_strvector(readcigar, &myreadcigar)) { return NULL; };
  SPLICING_FINALLY(splicing_strvector_destroy, &myreadcigar);
  
  SPLICING_PYCHECK(splicing_miso(mygff, gene, &myreadpos, 
				 (const char**) myreadcigar.table, 
				 readLength, noIterations, noBurnIn, noLag,
				 &myhyperp, &samples, &logLik, &rundata));

  splicing_vector_destroy(&myhyperp);
  splicing_vector_int_destroy(&myreadpos);
  splicing_strvector_destroy(&myreadcigar);
  SPLICING_FINALLY_CLEAN(3);
  
  result1=pysplicing_from_matrix(&samples);
  result2=pysplicing_from_vector(&logLik);
  result3=pysplicing_from_miso_rundata(&rundata);
  
  splicing_vector_destroy(&logLik);
  splicing_matrix_destroy(&samples);
  SPLICING_FINALLY_CLEAN(2);
  
  return Py_BuildValue("OOO", result1, result2, result3);
}

static PyObject* pysplicing_write_gff(PyObject *self, PyObject *args) {
  PyObject *gff;
  const char *filename;
  FILE *output;
  splicing_gff_t *mygff;
  
  if (!PyArg_ParseTuple(args, "Os", &gff, &filename)) { return NULL; }
  
  output = fopen(filename, "w");
  
  mygff = PyCObject_AsVoidPtr(gff);
  
  splicing_gff_write(output, mygff);
  
  fclose(output);
  
  Py_RETURN_NONE;
}

static PyObject* pysplicing_miso_paired(PyObject *self, PyObject*args) {
  PyObject *gff, *readpos, *readcigar, *hyperp=0;
  int gene, readLength, noIterations=5000, noBurnIn=500, noLag=10;
  double normalMean, normalVar, numDevs;
  splicing_gff_t *mygff;
  splicing_strvector_t myreadcigar;
  splicing_vector_int_t myreadpos;
  splicing_vector_t myhyperp;
  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_miso_rundata_t rundata;
  PyObject *result1, *result2, *result3;
  
  if (!PyArg_ParseTuple(args, "OiOOiddd|iiiO", &gff, &gene, &readpos, 
			&readcigar, &readLength, &normalMean, &normalVar,
			&numDevs, &noIterations, &noBurnIn, &noLag,
			&hyperp)) { return NULL; }

  mygff=PyCObject_AsVoidPtr(gff);
  
  SPLICING_PYCHECK(splicing_matrix_init(&samples, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &samples);
  SPLICING_PYCHECK(splicing_vector_init(&logLik, 0));
  SPLICING_FINALLY(splicing_vector_destroy, &logLik);
  
  if (pysplicing_to_vector_int(readpos, &myreadpos)) { return NULL; }
  SPLICING_FINALLY(splicing_vector_int_destroy, &myreadpos);
  
  if (hyperp) { 
    if (pysplicing_to_vector(hyperp, &myhyperp)) { return NULL; }
    SPLICING_FINALLY(splicing_vector_destroy, &myhyperp);
  } else {
    size_t i, noiso;
    SPLICING_PYCHECK(splicing_gff_noiso_one(mygff, gene, &noiso));
    SPLICING_PYCHECK(splicing_vector_init(&myhyperp, noiso));
    SPLICING_FINALLY(splicing_vector_destroy, &myhyperp);
    for (i=0; i<noiso; i++) { VECTOR(myhyperp)[i] = 1.0; }
  }
  if (pysplicing_to_strvector(readcigar, &myreadcigar)) { return NULL; }
  SPLICING_FINALLY(splicing_strvector_destroy, &myreadcigar);
  
  splicing_miso_paired(mygff, gene, &myreadpos,
		       (const char**) myreadcigar.table, readLength,
		       noIterations, noBurnIn, noLag, &myhyperp, 
		       /*insertProb=*/ 0, /*insertStart=*/ 0,
		       normalMean, normalVar, numDevs, &samples, &logLik,
		       &rundata);
  
  splicing_vector_destroy(&myhyperp);
  splicing_vector_int_destroy(&myreadpos);
  splicing_strvector_destroy(&myreadcigar);
  SPLICING_FINALLY_CLEAN(3);
  
  result1=pysplicing_from_matrix(&samples);
  result2=pysplicing_from_vector(&logLik);
  result3=pysplicing_from_miso_rundata(&rundata);
  
  splicing_vector_destroy(&logLik);
  splicing_matrix_destroy(&samples);
  SPLICING_FINALLY_CLEAN(2);
  
  return Py_BuildValue("OOO", result1, result2, result3);
}

static PyObject* pysplicing_create_gene(PyObject *self, PyObject *args) {
  PyObject *exons, *isoforms;
  const char *id="insilicogene", *seqid="seq1", *source="protein_coding";
  int strand=2;			/* unknown */
  splicing_vector_int_t myexons, myisoforms;
  splicing_gff_t *gff;
  
  if (!PyArg_ParseTuple(args, "OO|sssi", &exons, &isoforms, &id, &seqid, 
			&source, &strand)) { return NULL; }
  
  if (pysplicing_to_exons(exons, &myexons)) { return NULL; }
  SPLICING_FINALLY(splicing_vector_int_destroy, &myexons);
  if (pysplicing_to_isoforms(isoforms, &myisoforms)) { return NULL; }
  SPLICING_FINALLY(splicing_vector_int_destroy, &myisoforms);
  gff = malloc(sizeof(splicing_gff_t));
  if (!gff) { 
    splicing_error("Cannot create GFF", __FILE__, __LINE__, SPLICING_ENOMEM);
    splicingmodule_handle_splicing_error(); 
    return NULL; 
  } 
  SPLICING_FINALLY(splicing_free, gff);
  SPLICING_PYCHECK(splicing_gff_init(gff, 0));
  SPLICING_FINALLY(splicing_gff_destroy, gff);
    
  SPLICING_PYCHECK(splicing_create_gene(&myexons, &myisoforms, id, seqid, 
					source, strand, gff));

  splicing_vector_int_destroy(&myisoforms);
  splicing_vector_int_destroy(&myexons);
  SPLICING_FINALLY_CLEAN(4);
  
  return PyCObject_FromVoidPtr(gff, splicing_gff_destroy2);
}

static PyObject* pysplicing_simulate_reads(PyObject *self, PyObject *args) {
  PyObject *gff, *expression;
  int gene, noreads, readLength;
  splicing_vector_t myexpression;
  splicing_gff_t *mygff;
  splicing_vector_int_t isoform, position;
  splicing_strvector_t cigar;
  PyObject *risoform, *rposition, *rcigar;
  
  if (!PyArg_ParseTuple(args, "OiOii", &gff, &gene, &expression, &noreads,
			&readLength)) { return NULL; }

  mygff=PyCObject_AsVoidPtr(gff);
  if (pysplicing_to_vector(expression, &myexpression)) { return NULL; }
  SPLICING_FINALLY(splicing_vector_destroy, &myexpression);
  
  SPLICING_PYCHECK(splicing_strvector_init(&cigar, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &cigar);
  SPLICING_PYCHECK(splicing_vector_int_init(&position, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &position);
  SPLICING_PYCHECK(splicing_vector_int_init(&isoform, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isoform);
  
  SPLICING_PYCHECK(splicing_simulate_reads(mygff, gene, &myexpression, 
					   noreads, readLength,
					   &isoform, &position, &cigar));
  
  risoform = pysplicing_from_vector_int(&isoform);
  splicing_vector_int_destroy(&isoform);
  SPLICING_FINALLY_CLEAN(1);
  rposition = pysplicing_from_vector_int(&position);
  splicing_vector_int_destroy(&position);
  SPLICING_FINALLY_CLEAN(1);
  rcigar = pysplicing_from_strvector(&cigar);
  splicing_strvector_destroy(&cigar);
  SPLICING_FINALLY_CLEAN(1);
  
  splicing_vector_destroy(&myexpression);
  SPLICING_FINALLY_CLEAN(1);
  
  return Py_BuildValue("OOO", risoform, rposition, rcigar);
}

static PyObject* pysplicing_assignment_matrix(PyObject *self, 
					      PyObject *args) {
  PyObject *gff;
  int gene, readlength;
  splicing_matrix_t matrix;
  PyObject *rmatrix;
  splicing_gff_t *mygff;
  
  if (!PyArg_ParseTuple(args, "Oii", &gff, &gene, &readlength)) { 
    return NULL; 
  }
  
  mygff=PyCObject_AsVoidPtr(gff);
  SPLICING_PYCHECK(splicing_matrix_init(&matrix, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &matrix);
  
  SPLICING_PYCHECK(splicing_assignment_matrix(mygff, gene, readlength, 
					      &matrix));
  
  rmatrix = pysplicing_from_matrix(&matrix);
  splicing_matrix_destroy(&matrix);
  SPLICING_FINALLY_CLEAN(1);
  
  return Py_BuildValue("O", rmatrix);
}

static PyObject* pysplicing_gff_noiso(PyObject *self, PyObject *args) {
  
  PyObject *gff;
  splicing_gff_t *mygff;
  splicing_vector_int_t noiso;
  PyObject *rnoiso;
  
  if (!PyArg_ParseTuple(args, "O", &gff)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  SPLICING_PYCHECK(splicing_vector_int_init(&noiso, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &noiso);

  SPLICING_PYCHECK(splicing_gff_noiso(mygff, &noiso));
  
  rnoiso = pysplicing_from_vector_int(&noiso);
  splicing_vector_int_destroy(&noiso);
  SPLICING_FINALLY_CLEAN(1);

  return Py_BuildValue("O", rnoiso);
}  

static PyObject* pysplicing_gff_isolength(PyObject *self, PyObject *args) {

  PyObject *gff;
  splicing_gff_t *mygff;
  splicing_vector_int_t isolength, isolength_idx;
  PyObject *risolength;
  
  if (!PyArg_ParseTuple(args, "O", &gff)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  SPLICING_PYCHECK(splicing_vector_int_init(&isolength, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isolength);
  SPLICING_PYCHECK(splicing_vector_int_init(&isolength_idx, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isolength_idx);
  
  SPLICING_PYCHECK(splicing_gff_isolength(mygff, &isolength, &isolength_idx));
  
  risolength = pysplicing_from_vector_int_index(&isolength, &isolength_idx);

  splicing_vector_int_destroy(&isolength);
  splicing_vector_int_destroy(&isolength_idx);
  SPLICING_FINALLY_CLEAN(2);
  
  return Py_BuildValue("O", risolength);
}

static PyObject* pysplicing_solve_gene(PyObject *self, PyObject *args) {
  
  PyObject *gff, *readcigar, *position;
  int gene, readLength;
  splicing_gff_t *mygff;
  splicing_vector_int_t myposition;
  splicing_strvector_t myreadcigar;
  splicing_matrix_t match_matrix, assignment_matrix;
  splicing_vector_t expression;
  PyObject *r1, *r2, *r3;
  
  if (!PyArg_ParseTuple(args, "OiiOO", &gff, &gene, &readLength, &position,
			&readcigar)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  if (pysplicing_to_vector_int(position, &myposition)) { return NULL; }
  SPLICING_FINALLY(splicing_vector_int_destroy, &myposition);
  if (pysplicing_to_strvector(readcigar, &myreadcigar)) { return NULL; }
  SPLICING_FINALLY(splicing_strvector_destroy, &myreadcigar);

  SPLICING_PYCHECK(splicing_matrix_init(&match_matrix, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &match_matrix);
  SPLICING_PYCHECK(splicing_matrix_init(&assignment_matrix, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &assignment_matrix);
  SPLICING_PYCHECK(splicing_vector_init(&expression, 0));
  SPLICING_FINALLY(splicing_vector_destroy, &expression);

  SPLICING_PYCHECK(splicing_solve_gene(mygff, gene, readLength, &myposition, 
				       (const char **) myreadcigar.table,
				       &match_matrix, &assignment_matrix,
				       &expression));
  
  r1=pysplicing_from_matrix(&match_matrix);
  r2=pysplicing_from_matrix(&assignment_matrix);
  r3=pysplicing_from_vector(&expression);

  splicing_vector_destroy(&expression);
  splicing_matrix_destroy(&assignment_matrix);
  splicing_matrix_destroy(&match_matrix);
  splicing_strvector_destroy(&myreadcigar);
  splicing_vector_int_destroy(&myposition);
  SPLICING_FINALLY_CLEAN(5);

  return Py_BuildValue("OOO", r1, r2, r3);
}

static PyObject* pysplicing_simulate_paired_reads(PyObject *self,
						  PyObject *args) {
  PyObject *gff, *expression;
  int gene, noreads, readlength;
  double normalMean, normalVar, numDevs;
  splicing_gff_t *mygff;
  splicing_vector_t myexpression;
  splicing_vector_int_t isoform, position;
  splicing_strvector_t cigar;
  PyObject *r1, *r2, *r3;
  
  if (!PyArg_ParseTuple(args, "OiOiiddd", &gff, &gene, &expression,
			&noreads, &readlength, &normalMean, &normalVar,
			&numDevs)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  if (pysplicing_to_vector(expression, &myexpression)) { return NULL; }
  SPLICING_FINALLY(splicing_vector_destroy, &myexpression);
  
  SPLICING_PYCHECK(splicing_vector_int_init(&isoform, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &isoform);
  SPLICING_PYCHECK(splicing_vector_int_init(&position, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &position);
  SPLICING_PYCHECK(splicing_strvector_init(&cigar, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &cigar);
  
  SPLICING_PYCHECK(splicing_simulate_paired_reads(mygff, gene, &myexpression,
						  noreads, readlength,
						  /*insertProb=*/ 0,
						  /*insertStart=*/ 0,
						  normalMean, normalVar,
						  numDevs, &isoform, 
						  &position, &cigar));

  r1=pysplicing_from_vector_int(&isoform);
  r2=pysplicing_from_vector_int(&position);
  r3=pysplicing_from_strvector(&cigar);
  
  splicing_strvector_destroy(&cigar);
  splicing_vector_int_destroy(&position);
  splicing_vector_int_destroy(&isoform);
  splicing_vector_destroy(&myexpression);
  SPLICING_FINALLY_CLEAN(4);
  
  return Py_BuildValue("OOO", r1, r2, r3);
}

static PyObject* pysplicing_gene_complexity(PyObject *self, PyObject *args) {
  
  PyObject *gff;
  int gene, readlength, type, norm, paired;
  double normalMean, normalVar, numDevs;
  splicing_gff_t *mygff;
  double complexity;
  
  if (!PyArg_ParseTuple(args, "Oiiiiiddd", &gff, &gene, &readlength, 
			&type, &norm, &paired, &normalMean, &normalVar,
			&numDevs)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  
  SPLICING_PYCHECK(splicing_gene_complexity(mygff, gene, readlength,
					    type, norm, paired, 
					    /*insertProb=*/ 0, 0,
					    normalMean, normalVar, numDevs,
					    &complexity));
  
  return Py_BuildValue("d", complexity);
}
  
/* -------------------------------------------------------------------- */

static PyMethodDef pysplicing_methods[] = { 
  { "readGFF", pysplicing_read_gff, METH_VARARGS, "Read a GFF3 file." },
  { "writeGFF", pysplicing_write_gff, METH_VARARGS, "Write a GFF3 file." },
  { "createGene", pysplicing_create_gene, METH_VARARGS, "Create a gene." },
  { "simulateReads", pysplicing_simulate_reads, METH_VARARGS, 
    "Simulate reads from a gene." },
  { "assignmentMatrix", pysplicing_assignment_matrix, METH_VARARGS, 
    "Calculate the assignment matrix of a gene." },
  { "noIso", pysplicing_gff_noiso, METH_VARARGS, "Number if isoforms. " },
  { "isoLength", pysplicing_gff_isolength, METH_VARARGS, 
    "Length(s) of the various isoforms" },
  { "solveIsoGene", pysplicing_solve_gene, METH_VARARGS, 
    "Simple linear deconvolution of isoform expression for a single gene." },
  { "simulatePairedReads", pysplicing_simulate_paired_reads,
    METH_VARARGS, "Simulate paired end reads from a gene." },
  { "MISO"   , pysplicing_miso    , METH_VARARGS, "Run MISO." },
  { "MISOPaired", pysplicing_miso_paired, METH_VARARGS, 
    "MISO on paired-end data" },
  { "geneComplexity", pysplicing_gene_complexity, METH_VARARGS,
    "Gene complexity based on a linear model" },
  { NULL, NULL, 0, NULL }
};

extern PyObject* splicingmodule_InternalError;

PyMODINIT_FUNC initpysplicing(void) {
  PyObject *m = Py_InitModule3("pysplicing", pysplicing_methods,
			       "Python module for alternative splicing");
  if (m == NULL) { return; }
  Py_INCREF(m);

  splicingmodule_InternalError =
    PyErr_NewException("pysplicing.InternalError", PyExc_Exception, NULL);
  Py_INCREF(splicingmodule_InternalError);
  PyModule_AddObject(m, "InternalError", splicingmodule_InternalError);

  pysplicing_init_rng(m);
  splicing_set_error_handler(splicingmodule_splicing_error_hook);
  splicing_set_warning_handler(splicingmodule_splicing_warning_hook);
}
