##
## Sampling utilities
##
## Yarden Katz <yarden@mit.edu>
##
cimport numpy as np
import numpy as np

cimport cython

from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport log
from libc.math cimport sqrt
from libc.math cimport M_PI as PI
from libc.stdlib cimport RAND_MAX
from libc.stdlib cimport rand

cimport stat_helpers
cimport math_utils

cdef np.ndarray[double, ndim=1] \
  sample_indep_unit_normals(int N):
    """
    Draw a vector of N independent normal variables
    with unit variance (sigma^2 = 1).
    """
    cdef np.ndarray[double, ndim=1] samples = \
      np.empty(N, dtype=float)
    return samples

cdef double pyx_rand_normal_boxmuller():
    """
    Generate random sample from a unit normal, N(0, 1).
    Uses the non-polar form of Box-Muller transform.
    """
    cdef double U = 0.0
    cdef double V = 0.0
    cdef int phase = 0
    cdef double Z = 0.0

    if (phase == 0):
        U = (rand() + 1.) / (RAND_MAX + 2.)
        V = rand() / (RAND_MAX + 1.)
        Z = sqrt(-2 * log(U)) * sin(2 * PI * V)
    else:
        Z = sqrt(-2 * log(U)) * cos(2 * PI * V)
    phase = 1 - phase
    return Z

def py_rand_normal_boxmuller():
    return pyx_rand_normal_boxmuller()

    
