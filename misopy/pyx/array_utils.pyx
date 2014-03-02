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
