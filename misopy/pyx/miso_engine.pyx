##
## MISO engine
##
## Yarden Katz <yarden@mit.edu>
##
cimport numpy as np
import numpy as np

cimport matrix_utils
cimport array_utils
cimport math_utils
cimport sampling_utils


cdef class MISOEngine_PairedEnd:
  """
  MISO paired-end engine.
  """
  cdef int[:, :] reads

  
cdef class MISOEngine_SingleEnd:
  """
  MISO single-end engine.
  """
  cdef [:, :] reads
