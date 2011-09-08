
#include "random.h"
#include "splicing_random.h"
#include "splicing_error.h"

#include <limits.h>

/**
 * \ingroup python_interface_rng
 * \brief Internal data structure for storing references to the
 *        functions used from Python's random number generator.
 */
typedef struct {
  PyObject* randint_func;
  PyObject* random_func;
  PyObject* gauss_func;
} splicing_i_rng_Python_state_t;

static splicing_i_rng_Python_state_t splicing_rng_Python_state = {0, 0, 0};
static splicing_rng_t splicing_rng_Python = {0, 0, 0};

int splicing_rng_Python_init(void **state) {
  SPLICING_ERROR("Python RNG error, unsupported function called",
      SPLICING_EINTERNAL);
  return 0;
}

void splicing_rng_Python_destroy(void *state) {
  splicing_error("Python RNG error, unsupported function called",
      __FILE__, __LINE__, SPLICING_EINTERNAL);
}

/**
 * \ingroup python_interface_rng
 * \brief Sets the random number generator used by splicing.
 */
PyObject* splicing_rng_Python_set_generator(PyObject* self, PyObject* object) {
  splicing_i_rng_Python_state_t new_state, old_state;
  PyObject* func;

  if (object == Py_None) {
    /* Reverting to the default splicing random number generator instead
     * of the Python-based one */
    splicing_rng_set_default(&splicing_rng_default);
    Py_RETURN_NONE;
  }

#define GET_FUNC(name) {\
  func = PyObject_GetAttrString(object, name); \
  if (func == 0) \
    return NULL; \
  if (!PyCallable_Check(func)) {\
    PyErr_SetString(PyExc_TypeError, name "attribute must be callable"); \
    return NULL; \
  } \
}

  GET_FUNC("randint"); new_state.randint_func = func;
  GET_FUNC("random"); new_state.random_func = func;
  GET_FUNC("gauss"); new_state.gauss_func = func;

  old_state = splicing_rng_Python_state;
  splicing_rng_Python_state = new_state;
  Py_XDECREF(old_state.randint_func);
  Py_XDECREF(old_state.random_func);
  Py_XDECREF(old_state.gauss_func);

  splicing_rng_set_default(&splicing_rng_Python);

  Py_RETURN_NONE;
}

/**
 * \ingroup python_interface_rng
 * \brief Sets the seed of the random generator.
 */
int splicing_rng_Python_seed(void *state, unsigned long int seed) {
  SPLICING_ERROR("Python RNG error, unsupported function called",
      SPLICING_EINTERNAL);
  return 0;
}

/**
 * \ingroup python_interface_rng
 * \brief Generates an unsigned long integer using the Python random number generator.
 */
unsigned long int splicing_rng_Python_get(void *state) {
  PyObject* result = PyObject_CallFunction(splicing_rng_Python_state.randint_func, "kk", 0, LONG_MAX);
  unsigned long int retval;

  if (result == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    /* Fallback to the C random generator */
    return rand() * LONG_MAX;
  }
  retval = PyInt_AsLong(result);
  Py_DECREF(result);
  return retval;
}

/**
 * \ingroup python_interface_rng
 * \brief Generates a real number between 0 and 1 using the Python random number generator.
 */
double splicing_rng_Python_get_real(void *state) {
  PyObject* result = PyObject_CallFunction(splicing_rng_Python_state.random_func, NULL);
  double retval;

  if (result == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    /* Fallback to the C random generator */
    return rand();
  }

  retval = PyFloat_AsDouble(result);
  Py_DECREF(result);
  return retval;
}

/**
 * \ingroup python_interface_rng
 * \brief Generates a real number distributed according to the normal distribution
 *        around zero with unit variance.
 */
double splicing_rng_Python_get_norm(void *state) {
  PyObject* result = PyObject_CallFunction(splicing_rng_Python_state.gauss_func, "dd", 0.0, 1.0);
  double retval;

  if (result == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    /* Fallback to the C random generator */
    return 0;
  }

  retval = PyFloat_AsDouble(result);
  Py_DECREF(result);
  return retval;
}

/**
 * \ingroup python_interface_rng
 * \brief Specification table for Python's random number generator.
 *        This tells splicing which functions to call to obtain random numbers.
 */
splicing_rng_type_t splicing_rngtype_Python = {
  /* name= */      "Python random generator",
  /* min=  */      0,
  /* max=  */      LONG_MAX,
  /* init= */      splicing_rng_Python_init,
  /* destroy= */   splicing_rng_Python_destroy,
  /* seed= */      splicing_rng_Python_seed,
  /* get= */       splicing_rng_Python_get,
  /* get_real */   splicing_rng_Python_get_real,
  /* get_norm= */  splicing_rng_Python_get_norm,
  /* get_geom= */  0,
  /* get_binom= */ 0
};

void pysplicing_init_rng(PyObject* splicing_module) {
  PyObject* random_module;

  if (splicing_rng_Python.state != 0)
    return;

  random_module = PyImport_ImportModule("random");
  if (random_module == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    return;
  }

  splicing_rng_Python.type = &splicing_rngtype_Python;
  splicing_rng_Python.state = &splicing_rng_Python_state;

  if (splicing_rng_Python_set_generator(splicing_module, random_module) == 0) {
    PyErr_WriteUnraisable(PyErr_Occurred());
    PyErr_Clear();
    return;
  }
  Py_DECREF(random_module);
}
