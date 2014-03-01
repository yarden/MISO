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
cimport matrix_utils
cimport math_utils


##
## Generate samples from a unit normal distribution:
##   X ~ N(mu=0, sigma=1)
##
cdef double rand_normal_boxmuller():
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
    return rand_normal_boxmuller()


##
## Generate samples from N independent unit normal distributions
##
cdef np.ndarray[double, ndim=1] \
  sample_indep_unit_normals(int N):
    """
    Draw a row vector of N independent normal variables
    with mean = 0 and unit variance (sigma^2 = 1).
    """
    cdef np.ndarray[double, ndim=1] samples = \
      np.empty(N, dtype=float)
    cdef int i = 0
    for i in xrange(N):
        samples[i] = rand_normal_boxmuller()
    return samples


##
## Generate sample from a multivariate normal distribution
##
cdef np.ndarray[double, ndim=1] \
  sample_multivar_normal(np.ndarray[double, ndim=2] mu,
                         np.ndarray[double, ndim=2] L,
                         int k):
    """
    Draw a sample (vector) from a multivariate normal distribution
    given a vector mu and a matrix L which is the Cholesky
    decomposition of the covariance matrix Sigma of the
    distribution, i.e.
    
      N(mu, Sigma) where Sigma = LL'.

    If mu is a k-dimensional vector, the random sample from
    N(mu, Sigma) is generating by first drawing k-many random
    samples from independent unit normal distributions:

      Y = [y1, ..., yk]

    Then the random samples S are defined as:

      S = Ly + mu

    NOTE:
     (i) The argument L to this function is *not* the
    covariance matrix, but the Cholesky decomposition of it.

     (ii) mu must be a *column* vector! 
    """
    cdef int L_num_rows = L.shape[0]
    cdef int L_num_cols = L.shape[1]
    # Check that L is a square matrix
    if (L_num_rows != L_num_cols):
        print "Error: Cholesky decomposition matrix L must be square!"
    # Column vector S of samples to be generated and returned
    cdef np.ndarray[double, ndim=2] S = \
      np.empty((k, 1), dtype=float)
    # Column vector Y of independent random samples
    cdef np.ndarray[double, ndim=2] Y = \
      np.empty((k, 1), dtype=float)
    cdef int i = 0
    # Draw K-many independent samples from unit normal
    Y = matrix_utils.row_to_col_vect(sample_indep_unit_normals(k), k)
    # Generate the samples S = LY + mu
    # First compute S = LY
    print "MULTIPLYING: "
    print L
    print "with"
    print Y
    S = \
      matrix_utils.mat_times_mat(L, L_num_rows, L_num_cols,
                                 k, Y)
    # Now add mu: S = S + mu
    S = matrix_utils.mat_plus_mat(S, 2, k,
                                  mu, 2, k)
    return S


def py_sample_multivar_normal(np.ndarray[double, ndim=2] mu,
                              np.ndarray[double, ndim=2] L,
                              int k):
    return sample_multivar_normal(mu, L, k)


    
