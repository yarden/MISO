cimport cython
from cython.view cimport array as cvarray
from cpython.array cimport array, clone

cpdef array[double] get_double_array(int arr_size)
cpdef array[int] get_int_array(int arr_size)
