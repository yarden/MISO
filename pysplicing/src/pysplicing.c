
#include <Python.h>
#include <structmember.h>

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
  PyObject *gff, *reads, *hyperp=0;
  int gene, readLength, noIterations=5000, maxIterations=100000, 
    noBurnIn=500, noLag=10;
  int overhang=1;
  int no_chains=6;
  splicing_algorithm_t algo=SPLICING_ALGO_REASSIGN;
  splicing_miso_start_t start=SPLICING_MISO_START_AUTO;
  splicing_miso_stop_t stop=SPLICING_MISO_STOP_FIXEDNO;
  splicing_gff_t *mygff;
  splicing_replicate_reads_t myreads;
  splicing_miso_hyperprior_t hyperprior;
  double logistic_mean = 0.0, logistic_var = 3.0;
  double replicate_mean_prior_mean = 0.0;
  double replicate_mean_prior_var = 100.0;
  double replicate_var_prior_numobs = 1/10.0;
  double replicate_var_prior_var = 1/10.0;
  splicing_matrix_t pop_samples;
  splicing_vector_ptr_t samples;	/* matrix_t */
  splicing_vector_t logLik;
  splicing_matrix_t class_templates;
  splicing_vector_ptr_t class_counts; /* vector_t */
  splicing_vector_ptr_t assignment;   /* vector_int_t */
  splicing_miso_rundata_t rundata;
  PyObject *r1, *r2, *r3, *r4, *r5, *r6;

  hyperprior.prior = SPLICING_MISO_PRIOR_DIRICHLET;

  if (!PyArg_ParseTuple(args,
			"O"	/* gff */
			"i"	/* gene */
			"O"	/* reads */
			"i"	/* readLength */
			"|i"	/* noIterations */
			"i"	/* noBurnIn */
			"i"	/* noLag */
			"i" 	/* prior */
			"O"	/* hyperp */
			"d"	/* logistic_mean */
			"d"	/* logistic_var */
			"d"     /* replicate_mean_prior_mean */
			"d"     /* replicate_mean_prior_var */
			"d"     /* replicate_var_prior_numobs */
			"d" 	/* replicate_var_prior_var */
			"i"	/* overhang */
			"i"	/* no_chains */
			"i"	/* start */
			"i"	/* stop */
			"i",	/* algo */
			&gff, &gene, &reads,
			&readLength, &noIterations, &noBurnIn, &noLag,
			&hyperprior.prior, &hyperp,
			&logistic_mean, &logistic_var,
			&replicate_mean_prior_mean,
			&replicate_mean_prior_var,
			&replicate_var_prior_numobs,
			&replicate_var_prior_var,
			&overhang, &no_chains, &start, &stop, &algo)) {
    return NULL; 
  }
  
  mygff=PyCObject_AsVoidPtr(gff);

  SPLICING_PYCHECK(splicing_matrix_init(&pop_samples, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &pop_samples);
  SPLICING_PYCHECK(splicing_vector_ptr_init(&samples, 0));
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &samples);
  SPLICING_PYCHECK(splicing_vector_init(&logLik, 0));
  SPLICING_FINALLY(splicing_vector_destroy, &logLik);
  SPLICING_PYCHECK(splicing_matrix_init(&class_templates, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &class_templates);
  SPLICING_PYCHECK(splicing_vector_ptr_init(&class_counts, 0));
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &class_counts);
  SPLICING_PYCHECK(splicing_vector_ptr_init(&assignment, 0));
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &assignment);
  if (pysplicing_to_replicate_reads(reads, &myreads)) { return NULL; }
  SPLICING_FINALLY(splicing_replicate_reads_destroy, &myreads);

  if (hyperp) { 
    if (pysplicing_to_vector(hyperp, &hyperprior.dirichlet_hyperp)) { return NULL; }
    SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.dirichlet_hyperp);
  } else {
    size_t i, noiso;
    SPLICING_PYCHECK(splicing_gff_noiso_one(mygff, gene, &noiso));
    SPLICING_PYCHECK(splicing_vector_init(&hyperprior.dirichlet_hyperp, noiso));
    SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.dirichlet_hyperp);
    for (i=0; i<noiso; i++) { VECTOR(hyperprior.dirichlet_hyperp)[i] = 1.0; }
  }

  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_mean, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_mean);
  splicing_vector_fill(&hyperprior.logistic_mean, logistic_mean);
  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_var, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_var);
  splicing_vector_fill(&hyperprior.logistic_var, logistic_var);

  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_mean_mean, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_mean_mean);
  splicing_vector_fill(&hyperprior.logistic_mean_mean, replicate_mean_prior_mean);
  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_mean_var, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_mean_var);
  splicing_vector_fill(&hyperprior.logistic_mean_var, replicate_mean_prior_var);

  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_var_numobs, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_var_numobs);
  splicing_vector_fill(&hyperprior.logistic_var_numobs, replicate_var_prior_numobs);
  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_var_var, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_var_var);
  splicing_vector_fill(&hyperprior.logistic_var_var, replicate_var_prior_var);
  
  SPLICING_PYCHECK(splicing_miso(mygff, gene, &myreads,
				 readLength, overhang, no_chains,
				 noIterations, maxIterations, 
				 noBurnIn, noLag, &hyperprior, algo,
				 start, stop, 0, &pop_samples, &samples, &logLik,
				 /*match_matrix=*/ 0, &class_templates,
				 &class_counts, &assignment, &rundata));

  splicing_vector_destroy(&hyperprior.logistic_var_var);
  splicing_vector_destroy(&hyperprior.logistic_var_numobs);
  splicing_vector_destroy(&hyperprior.logistic_mean_var);
  splicing_vector_destroy(&hyperprior.logistic_mean_mean);
  splicing_vector_destroy(&hyperprior.logistic_mean);
  splicing_vector_destroy(&hyperprior.logistic_var);
  SPLICING_FINALLY_CLEAN(6);

  splicing_vector_destroy(&hyperprior.dirichlet_hyperp);
  splicing_replicate_reads_destroy(&myreads);
  SPLICING_FINALLY_CLEAN(2);
  
  r6=pysplicing_from_miso_rundata(&rundata);

  r5=pysplicing_from_vectorlist_int(&assignment);
  splicing_vector_ptr_destroy(&assignment); SPLICING_FINALLY_CLEAN(1);

  r4=pysplicing_from_vectorlist(&class_counts);
  splicing_vector_ptr_destroy(&class_counts); SPLICING_FINALLY_CLEAN(1);

  splicing_matrix_transpose(&class_templates);
  r3=pysplicing_from_matrix(&class_templates);
  splicing_matrix_destroy(&class_templates); SPLICING_FINALLY_CLEAN(1);

  r2=pysplicing_from_vector(&logLik);
  splicing_vector_destroy(&logLik); SPLICING_FINALLY_CLEAN(1);

  r1=pysplicing_from_matrix(&pop_samples);
  splicing_matrix_destroy(&pop_samples); SPLICING_FINALLY_CLEAN(1);

  splicing_vector_ptr_destroy(&samples); SPLICING_FINALLY_CLEAN(1);
  
  return Py_BuildValue("OOOOOO", r1, r2, r3, r4, r5, r6);
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
  PyObject *gff, *reads, *hyperp=0;
  int gene, readLength, noIterations=5000, maxIterations=100000, 
    noBurnIn=500, noLag=10;
  int overhang=1;
  int no_chains=6;
  splicing_miso_start_t start=SPLICING_MISO_START_AUTO;
  splicing_miso_stop_t stop=SPLICING_MISO_STOP_FIXEDNO;
  double normalMean, normalVar, numDevs;
  splicing_gff_t *mygff;
  splicing_replicate_reads_t myreads;
  splicing_miso_hyperprior_t hyperprior;
  double logistic_mean = 0.0, logistic_var = 3.0;
  double replicate_mean_prior_mean = 0.0;
  double replicate_mean_prior_var = 100.0;
  double replicate_var_prior_numobs = 1/10.0;
  double replicate_var_prior_var = 1/10.0;
  splicing_matrix_t pop_samples;
  splicing_vector_ptr_t samples; /* matrix_t */
  splicing_vector_t logLik;
  splicing_matrix_t bin_class_templates;
  splicing_vector_ptr_t bin_class_counts;
  splicing_vector_ptr_t assignment;
  splicing_miso_rundata_t rundata;
  PyObject *r1, *r2, *r3, *r4, *r5, *r6;

  hyperprior.prior = SPLICING_MISO_PRIOR_DIRICHLET;

  if (!PyArg_ParseTuple(args,
			"O"	/* gff */
			"i"	/* gene */
			"O"	/* reads */
			"i"	/* readLength */
			"d"	/* normalMean */
			"d"	/* normalVar */
			"d"	/* numDevs */
			"|i"	/* noIterations */
			"i"	/* noBurnIn */
			"i"	/* noLag */
			"i"     /* prior */
			"O"	/* hyperp */
			"d"	/* logistic_mean */
			"d"	/* logistic_var */
			"d"     /* replicate_mean_prior_mean */
			"d"     /* replicate_mean_prior_var */
			"d"     /* replicate_var_prior_numobs */
			"d" 	/* replicate_var_prior_var */
			"i"	/* overhang */
			"i"	/* no_chains */
			"i"	/* start */
			"i",	/* stop */
			&gff, &gene, &reads,
			&readLength, &normalMean, &normalVar,
			&numDevs, &noIterations, &noBurnIn, &noLag,
			&hyperprior.dirichlet_hyperp, &hyperp,
			&logistic_mean, &logistic_var,
			&replicate_mean_prior_mean,
			&replicate_mean_prior_var,
			&replicate_var_prior_numobs,
			&replicate_var_prior_var,
			&overhang, &no_chains, &start, &stop)) {
    return NULL; 
  }

  mygff=PyCObject_AsVoidPtr(gff);
  
  SPLICING_PYCHECK(splicing_matrix_init(&pop_samples, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &pop_samples);
  SPLICING_PYCHECK(splicing_vector_ptr_init(&samples, 0));
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &samples);
  SPLICING_PYCHECK(splicing_vector_init(&logLik, 0));
  SPLICING_FINALLY(splicing_vector_destroy, &logLik);
  SPLICING_PYCHECK(splicing_matrix_init(&bin_class_templates, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &bin_class_templates);
  SPLICING_PYCHECK(splicing_vector_ptr_init(&bin_class_counts, 0));
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &bin_class_counts);
  SPLICING_PYCHECK(splicing_vector_ptr_init(&assignment, 0));
  SPLICING_FINALLY(splicing_vector_ptr_destroy, &assignment);
  
  if (pysplicing_to_replicate_reads(reads, &myreads)) { return NULL; }
  SPLICING_FINALLY(splicing_replicate_reads_destroy, &myreads);
  
  if (hyperp) { 
    if (pysplicing_to_vector(hyperp, &hyperprior.dirichlet_hyperp)) { return NULL; }
    SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.dirichlet_hyperp);
  } else {
    size_t i, noiso;
    SPLICING_PYCHECK(splicing_gff_noiso_one(mygff, gene, &noiso));
    SPLICING_PYCHECK(splicing_vector_init(&hyperprior.dirichlet_hyperp, noiso));
    SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.dirichlet_hyperp);
    for (i=0; i<noiso; i++) { VECTOR(hyperprior.dirichlet_hyperp)[i] = 1.0; }
  }

  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_mean, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_mean);
  splicing_vector_fill(&hyperprior.logistic_mean, logistic_mean);
  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_var, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_var);
  splicing_vector_fill(&hyperprior.logistic_var, logistic_var);

  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_mean_mean, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_mean_mean);
  splicing_vector_fill(&hyperprior.logistic_mean_mean, replicate_mean_prior_mean);
  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_mean_var, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_mean_var);
  splicing_vector_fill(&hyperprior.logistic_mean_var, replicate_mean_prior_var);

  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_var_numobs, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_var_numobs);
  splicing_vector_fill(&hyperprior.logistic_var_numobs, replicate_var_prior_numobs);
  SPLICING_PYCHECK(splicing_vector_init(&hyperprior.logistic_var_var, no_chains));
  SPLICING_FINALLY(splicing_vector_destroy, &hyperprior.logistic_var_var);
  splicing_vector_fill(&hyperprior.logistic_var_var, replicate_var_prior_var);
  
  splicing_miso_paired(mygff, gene, &myreads,
		       readLength,
		       overhang, no_chains, noIterations, maxIterations,
		       noBurnIn, noLag, &hyperprior,
		       start, stop, /*start_psi=*/ 0,
		       /*insertProb=*/ 0, /*insertStart=*/ 0,
		       normalMean, normalVar, numDevs, &pop_samples,
		       &samples, &logLik,
		       /*match_matrix=*/ 0, /*class_templates=*/ 0, 
		       /*class_counts=*/ 0, &bin_class_templates, 
		       &bin_class_counts, &assignment, &rundata);

  splicing_vector_destroy(&hyperprior.logistic_var_var);
  splicing_vector_destroy(&hyperprior.logistic_var_numobs);
  splicing_vector_destroy(&hyperprior.logistic_mean_var);
  splicing_vector_destroy(&hyperprior.logistic_mean_mean);
  splicing_vector_destroy(&hyperprior.logistic_mean);
  splicing_vector_destroy(&hyperprior.logistic_var);
  SPLICING_FINALLY_CLEAN(6);
  
  splicing_vector_destroy(&hyperprior.dirichlet_hyperp);
  splicing_replicate_reads_destroy(&myreads);
  SPLICING_FINALLY_CLEAN(2);
  
  r6=pysplicing_from_miso_rundata(&rundata);

  r5=pysplicing_from_vectorlist_int(&assignment);
  splicing_vector_ptr_destroy(&assignment); SPLICING_FINALLY_CLEAN(1);

  r4=pysplicing_from_vectorlist(&bin_class_counts);
  splicing_vector_ptr_destroy(&bin_class_counts); SPLICING_FINALLY_CLEAN(1);

  splicing_matrix_transpose(&bin_class_templates);
  r3=pysplicing_from_matrix(&bin_class_templates);
  splicing_matrix_destroy(&bin_class_templates); SPLICING_FINALLY_CLEAN(1);

  r2=pysplicing_from_vector(&logLik);
  splicing_vector_destroy(&logLik); SPLICING_FINALLY_CLEAN(1);

  r1=pysplicing_from_matrix(&pop_samples);
  splicing_matrix_destroy(&pop_samples); SPLICING_FINALLY_CLEAN(1);

  splicing_vector_ptr_destroy(&samples);
  
  return Py_BuildValue("OOOOOO", r1, r2, r3, r4, r5, r6);
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
					   &isoform, &position, &cigar, 
					   /*sample_prob=*/ 0));
  
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
  int gene, readlength, overhang=1;
  splicing_matrix_t matrix;
  PyObject *rmatrix;
  splicing_gff_t *mygff;
  
  if (!PyArg_ParseTuple(args, "Oii|i", &gff, &gene, &readlength, 
			&overhang)) { 
    return NULL; 
  }
  
  mygff=PyCObject_AsVoidPtr(gff);
  SPLICING_PYCHECK(splicing_matrix_init(&matrix, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &matrix);
  
  SPLICING_PYCHECK(splicing_assignment_matrix(mygff, gene, readlength, 
					      overhang, &matrix));
  
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
  int gene, readLength, overhang=1;
  splicing_gff_t *mygff;
  splicing_vector_int_t myposition;
  splicing_strvector_t myreadcigar;
  splicing_matrix_t match_matrix, assignment_matrix;
  splicing_vector_t expression;
  PyObject *r1, *r2, *r3;
  
  if (!PyArg_ParseTuple(args, "OiiOO|i", &gff, &gene, &readLength, &position,
			&readcigar, &overhang)) { return NULL; }
  
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

  SPLICING_PYCHECK(splicing_solve_gene(mygff, gene, readLength, overhang,
				       &myposition, 
				       (const char **) myreadcigar.table,
				       &match_matrix, /*nomatch=*/ 0, 
				       &assignment_matrix,
				       &expression, /*residuals=*/ 0, 
				       /*scale=*/ 0));
  
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
						  &position, &cigar, 0));

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
  int gene, readlength, type, norm, paired, overhang=1;
  double normalMean, normalVar, numDevs;
  splicing_gff_t *mygff;
  double complexity;
  
  if (!PyArg_ParseTuple(args, "Oiiiiiddd|i", &gff, &gene, &readlength, 
			&type, &norm, &paired, &normalMean, &normalVar,
			&numDevs, &overhang)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  
  SPLICING_PYCHECK(splicing_gene_complexity(mygff, gene, readlength,
					    overhang, type, norm, paired, 
					    /*insertProb=*/ 0, 0,
					    normalMean, normalVar, numDevs,
					    &complexity));
  
  return Py_BuildValue("d", complexity);
}

static PyObject* pysplicing_gff_nogenes(PyObject *self, PyObject *args) {
  
  PyObject *gff;
  size_t nogenes;
  int nogenes2;
  splicing_gff_t *mygff;
  
  if (!PyArg_ParseTuple(args, "O", &gff)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  
  SPLICING_PYCHECK(splicing_gff_nogenes(mygff, &nogenes));
  nogenes2 = (int) nogenes;		/* in case int and size_t are different */
    
  return Py_BuildValue("i", nogenes2);
} 

static PyObject* pysplicing_from_gff(PyObject *self, PyObject *args) {
  
  PyObject *gff;
  splicing_gff_t *mygff;
  PyObject *seqids, *sources, *genes, *transcripts, *seqid, *source,
    *type, *start, *end, *score, *strand, *phase, *ID, *parent;
  
  if (!PyArg_ParseTuple(args, "O", &gff)) { return NULL; }
  
  mygff=PyCObject_AsVoidPtr(gff);
  
  seqids=pysplicing_from_strvector(&mygff->seqids);
  sources=pysplicing_from_strvector(&mygff->sources);
  genes=pysplicing_from_vector_int(&mygff->genes);
  transcripts=pysplicing_from_vector_int(&mygff->transcripts);
  seqid=pysplicing_from_vector_int(&mygff->seqid);
  source=pysplicing_from_vector_int(&mygff->source);
  type=pysplicing_from_vector_int(&mygff->type);
  start=pysplicing_from_vector_int(&mygff->start);
  end=pysplicing_from_vector_int(&mygff->end);
  score=pysplicing_from_vector(&mygff->score);
  strand=pysplicing_from_vector_int(&mygff->strand);
  phase=pysplicing_from_vector_int(&mygff->phase);
  ID=pysplicing_from_strvector(&mygff->ID);
  parent=pysplicing_from_vector_int(&mygff->parent);
 
  return Py_BuildValue("OOOOOOOOOOOOOO", seqids, sources, genes, transcripts,
		       seqid, source, type, start, end, score, strand, 
		       phase, ID, parent);
}

/* TODO: Check that it works for missing values */

static PyObject* pysplicing_to_gff(PyObject *self, PyObject *args) {
  
  PyObject *pygff, *entries;
  size_t i, noRec;
  splicing_gff_t *cgff;
  PyObject *IDkey, *Parentkey;
  
  if (!PyArg_ParseTuple(args, "O", &pygff)) { return NULL; }

  if (!PyObject_HasAttrString(pygff, "_GFFDatabase__entries")) {
    splicingmodule_handle_splicing_error();
    return NULL;
  }
  entries=PyObject_GetAttrString(pygff, "_GFFDatabase__entries");
  noRec=PySequence_Size(entries);

  IDkey=PyString_FromString("ID");
  Parentkey=PyString_FromString("Parent");

  cgff=malloc(sizeof(splicing_gff_t));
  if (!cgff) { splicingmodule_handle_splicing_error(); return NULL; }
  splicing_gff_init(cgff, noRec);
  
  for (i=0; i<noRec; i++) {
    PyObject *rec=0, *seqid=0, *source=0, *type=0, *start=0, *end=0,
      *score=0, *strand=0, *phase=0, *attributes=0, *ID=0, *Parent=0;
    char *Cseqid, *Csource, *CID, *Cparent=0, *Ctype;
    splicing_type_t Ctype2;
    int Cstart, Cend, Cphase;
    double Cscore;
    splicing_strand_t Cstrand;
    
    rec=PySequence_GetItem(entries, i);
    seqid=PyObject_GetAttrString(rec, "seqid");
    source=PyObject_GetAttrString(rec, "source");
    type=PyObject_GetAttrString(rec, "type");
    start=PyObject_GetAttrString(rec, "start");
    end=PyObject_GetAttrString(rec, "end");
    score=PyObject_GetAttrString(rec, "score");
    strand=PyObject_GetAttrString(rec, "strand");
    phase=PyObject_GetAttrString(rec, "phase");
    attributes=PyObject_GetAttrString(rec, "attributes");
    ID=PyDict_GetItem(attributes, IDkey);
    Parent=PyDict_GetItem(attributes, Parentkey);

    Cseqid=PyString_AsString(seqid);
    Csource=PyString_AsString(source);

    Ctype=PyString_AsString(type);
    if (!strcmp(Ctype, "gene")) { 
      Ctype2=SPLICING_TYPE_GENE;
    } else if (!strcmp(Ctype, "mRNA")) {
      Ctype2=SPLICING_TYPE_MRNA;
    } else if (!strcmp(Ctype, "exon")) {
      Ctype2=SPLICING_TYPE_EXON;
    } else if (!strcmp(Ctype, "CDS")) {
      Ctype2=SPLICING_TYPE_CDS;
    } else if (!strcmp(Ctype, "start_codon")) {
      Ctype2=SPLICING_TYPE_START_CODON;
    } else if (!strcmp(Ctype, "stop_codon")) {
      Ctype2=SPLICING_TYPE_STOP_CODON;
    } else {
      splicing_error("Cannot create GFF", __FILE__, __LINE__,
		     SPLICING_EINVAL);
      splicingmodule_handle_splicing_error();
      return NULL;
    }

    Cstart = (int) PyInt_AsLong(start);
    Cend = (int) PyInt_AsLong(end);
    Cscore=PyFloat_AsDouble(score);
    Cstrand = (int) PyInt_AsLong(strand);
    Cphase = (int) PyInt_AsLong(phase);
    CID=PyString_AsString(ID);
    if (Parent) { Cparent=PyString_AsString(Parent); }

    SPLICING_PYCHECK(splicing_gff_append(cgff, Cseqid, Csource, Ctype2, 
					 Cstart, Cend, Cscore, Cstrand,
					 Cphase, CID, Cparent));

    Py_DECREF(rec); Py_DECREF(seqid); Py_DECREF(source); Py_DECREF(type);
    Py_DECREF(start); Py_DECREF(end); Py_DECREF(score); Py_DECREF(strand);
    Py_DECREF(phase); Py_DECREF(attributes);
  }

  Py_DECREF(entries);
  Py_DECREF(IDkey);
  Py_DECREF(Parentkey);

  SPLICING_PYCHECK(splicing_gff_reindex(cgff));
  
  return PyCObject_FromVoidPtr(cgff, splicing_gff_destroy2);
}
  
/* -------------------------------------------------------------------- */
/* These are internal, and only used for testing the C functions.       */

static PyObject* pysplicing_test_rng_get_invchi2(PyObject *self, PyObject *args) {

  splicing_vector_t nums;
  int n = 1000;
  double median, iqr;

  SPLICING_PYCHECK(splicing_vector_init(&nums, 1000));

  /* Generate 1000 numbers, and take their meadian and iqr */
#define TEST(nu, tau2) {					\
    int i;							\
    for (i = 0; i < n; i++) {					\
      VECTOR(nums)[i] = RNG_INVCHI2((nu), (tau2));		\
    }								\
    median = splicing_vector_median(&nums);			\
    iqr    = splicing_vector_iqr(&nums);			\
    /* printf("%f %f %f %f\n", nu, tau2, median, iqr);	*/	\
  }

  TEST(1.0, 0.5);
  if (median < 0.8 || median > 1.5) return Py_False;
  if (iqr    < 3.0 || iqr    > 7.0) return Py_False;

  TEST(1.0, 1.0);
  if (median < 1.5 || median >  3.0) return Py_False;
  if (iqr    < 5.0 || iqr    > 20.0) return Py_False;

  TEST(2.0, 0.5);
  if (median < 0.6 || median > 0.9) return Py_False;
  if (iqr    < 1.0 || iqr    > 2.0) return Py_False;

  TEST(2.0, 1.0);
  if (median < 1.2 || median > 1.8) return Py_False;
  if (iqr    < 2.0 || iqr    > 3.7) return Py_False;

  splicing_vector_destroy(&nums);
  return Py_True;
}

static PyObject* pysplicing_test_vector_median(PyObject *self, PyObject *args) {

  PyObject *v;
  splicing_vector_t myv;
  double median;

  if (!PyArg_ParseTuple(args, "O", &v)) { return NULL; }
  SPLICING_PYCHECK(pysplicing_to_vector(v, &myv));

  median = splicing_vector_median(&myv);

  return PyFloat_FromDouble(median);
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
  { "noGenes", pysplicing_gff_nogenes, METH_VARARGS, 
    "Number of genes in a GFF object." },
  { "i_fromGFF", pysplicing_from_gff, METH_VARARGS, 
    "Convert a C GFF structure to Python" },
  { "toGFF", pysplicing_to_gff, METH_VARARGS, "Convert a Python GFF to C" },

  /* These are internal and used for testing the C functions */

  { "i_test_rng_get_invchi2", pysplicing_test_rng_get_invchi2, METH_VARARGS,
    "Test drawing from a scaled inverse Chi squared distribution" },
  { "i_test_vector_median", pysplicing_test_vector_median, METH_VARARGS,
    "Test median of a vector" },

  { NULL, NULL, 0, NULL }
};

extern PyObject* splicingmodule_InternalError;

PyMODINIT_FUNC initpysplicing(void) {

  PyObject *m = Py_InitModule3("pysplicing", pysplicing_methods,
			       "Python module for alternative splicing");

  if (m == NULL) { return; }
  Py_INCREF(m);

  /* New exceptions */

  splicingmodule_InternalError =
    PyErr_NewException("pysplicing.InternalError", PyExc_Exception, NULL);
  Py_INCREF(splicingmodule_InternalError);
  PyModule_AddObject(m, "InternalError", splicingmodule_InternalError);

  /* Some callbacks */

  pysplicing_init_rng(m);
  splicing_set_error_handler(splicingmodule_splicing_error_hook);
  splicing_set_warning_handler(splicingmodule_splicing_warning_hook);

}
