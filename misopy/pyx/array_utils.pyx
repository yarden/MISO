##
## Utilities for working with Cython arrays
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython
from cython.view cimport array as cvarray
from cpython.array cimport array, clone


cpdef array[double] get_double_array(int arr_size):
    """
    Return a clone of a double array. Used for fast
    creation of Cython arrays.
    """
    DOUBLE_ARRAY_1D = array("d")
    arr = clone(DOUBLE_ARRAY_1D, arr_size, False)
    return arr


cpdef array[int] get_int_array(int arr_size):
    """
    Return a clone of an int array. Used for fast creation
    of Cython arrays.
    """
    INT_ARRAY_1D = array("i")
    arr = clone(INT_ARRAY_1D, arr_size, False)
    return arr


cpdef double \
  sum_array(double[:] input_array,
            int array_len):
    """
    Compute sum of elements of doubles array.
    
    Args:
      input_array: double array (1d)

      array_len: Length of array
      
    Returns:
      Sum of array
    """
    cdef int j = 0
    cdef double result = 0.0
    for j in xrange(array_len):
        result += input_array[j]
    return result


