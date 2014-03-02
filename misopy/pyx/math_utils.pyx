##
## Math utilities
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython

from libc.math cimport log
from libc.math cimport exp

cimport array_utils


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef double my_logsumexp(double[:] log_vector,
                         int vector_len):
    """
    Log sum exp.

    Parameters:
    -----------

    log_vector : array of floats corresponding to log values.
    vector_len : int, length of vector.

    Returns:
    --------

    Result of log(sum(exp(log_vector)))
    """
    cdef double curr_exp_value = 0.0
    cdef double sum_of_exps = 0.0
    cdef double log_sum_of_exps = 0.0
    cdef int curr_elt = 0
    # First find the maximum value first
    cdef double max_val = log_vector[0]
    for curr_elt in xrange(vector_len):
        if (log_vector[curr_elt] > max_val):
            max_val = log_vector[curr_elt]
    # Subtract maximum value from the rest
    for curr_elt in xrange(vector_len):
        curr_exp_value = exp(log_vector[curr_elt] - max_val)
        sum_of_exps += curr_exp_value
    # Now take log of the sum of exp values and add
    # back the missing value
    log_sum_of_exps = log(sum_of_exps) + max_val
    return log_sum_of_exps


cpdef double[:] logit(double[:] p,
                      int p_len):
    """
    Logit transform.
    
    Takes a vector with values \in (0, 1) and transforms each \
    to values \in (-inf, inf).
    """
    cdef double[:] logit_vals = array_utils.get_double_array(p_len)
    cdef int i = 0
    for i in xrange(p_len):
        logit_vals[i] = log(p[i] / (1-p[i]))
    return logit_vals


cpdef double[:] logit_inv(double[:] p,
                          int p_len):
    """
    Inverse Logit transform (logit-1()).
    """
    # original python code
    #denom = append(x, 0)
    #p = exp(x)/(sum(exp(denom)))
    cdef double[:] logit_inv_vals = array_utils.get_double_array(p_len)
    # The denominator is sum(exp(x)) + exp(0)
    cdef double denom = exp(0)
    cdef int i = 0
    # Calculate denominator first
    for i in xrange(p_len):
        denom += exp(p[i])
    # Now calculate inverse logit using denominator
    for i in xrange(p_len):
        logit_inv_vals[i] = exp(p[i]) / denom
    return logit_inv_vals

