##
## MISO engine
##
cimport numpy as np
import numpy as np

import misopy
import misopy.matrix_utils as matrix_utils

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_float_t


cdef class MISOEngine_PairedEnd:
  """
  MISO paired-end engine.
  """
  cdef np.ndarray[DTYPE_t, ndim=2] reads

  
cdef class MISOEngine_SingleEnd:
  """
  MISO single-end engine.
  """
  cdef np.ndarray[DTYPE_t, ndim=2] reads
