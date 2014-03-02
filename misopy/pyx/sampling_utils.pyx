##
## Sampling utilities
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython

from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport log
from libc.math cimport sqrt
from libc.math cimport M_PI as PI
from libc.stdlib cimport RAND_MAX
from libc.stdlib cimport rand

from cython.view cimport array as cvarray

cimport stat_helpers
cimport array_utils
cimport matrix_utils
cimport math_utils


##
## Generate samples from a unit normal distribution:
##   X ~ N(mu=0, sigma=1)
##
cpdef double rand_normal_boxmuller():
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


##
## Generate samples from N independent unit normal distributions
##
cdef double[:] \
  sample_indep_unit_normals(int N):
    """
    Draw a row vector of N independent normal variables
    with mean = 0 and unit variance (sigma^2 = 1).
    """
    cdef double[:] samples = \
      array_utils.get_double_array(N)
    cdef int i = 0
    for i in xrange(N):
        samples[i] = rand_normal_boxmuller()
    return samples


##
## Generate sample from a multivariate normal distribution
##
cpdef double[:] \
  sample_multivar_normal(double[:] mu,
                         double[:, :] L,
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

     (iii) The function returns a *row* vector, S'
    """
    cdef int L_num_rows = L.shape[0]
    cdef int L_num_cols = L.shape[1]
    # Check that L is a square matrix
    if (L_num_rows != L_num_cols):
        print "Error: Cholesky decomposition matrix L must be square!"
    cdef double[:] S_prod = \
      array_utils.get_double_array(k)
    # Column vector S_prod of samples to be generated and returned
    cdef double[:] S_final = \
      array_utils.get_double_array(k)
    # Column vector Y of independent random samples
    cdef double[:] Y = \
      array_utils.get_double_array(k)
    cdef int i = 0
    # Draw K-many independent samples from unit normal
    Y = sample_indep_unit_normals(k)
    # Generate the samples S = LY + mu
    # First compute S = LY, get result as 1d array
    S_prod = \
      matrix_utils.mat_times_col_vect(L,
                                      L_num_rows,
                                      L_num_cols,
                                      Y,
                                      k,
                                      S_prod)
    # Now add mu: LY + mu
    S_final = matrix_utils.vect_plus_vect(S_prod, mu, k, S_final)
    return S_final

    



    
