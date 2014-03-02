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


cdef class MISO_SingleEnd:
  """
  MISO single-end engine.
  """
  ## Gene information
  # Number of isoforms
  cdef int num_iso
  # Isoform lengths
  cdef int[:] iso_lens 
  
  ## Reads information
  cdef int[:, :] reads
  cdef int read_len
  cdef int overhang_len

  # Sampler parameters
  cdef int num_iters
  cdef int burn_in
  cdef int lag
  # Drift proposal distribution's Sigma value
  cdef double proposal_sigma
  # Covariance matrix of proposal distribution
  cdef double[:, :] covar_mat
  # Cholesky decomposition of covar matrix
  cdef double[:, :] covar_L
  
  def __cinit__(self):
      """
      Initialize variables prior to Python.
      """
      pass

  
  def __init__(self, reads, iso_lens, read_len,
               covar_mat,
               covar_L,
               overhang_len=1,
               num_iters=2500,
               burn_in=500
               lag=10):
      """
      Initialize parameters for Python instantiation of object.
      """
      self.reads = reads
      self.read_len = read_len
      self.num_iso = iso_lens.shape[0]
      self.covar_mat = covar_mat
      self.covar_L = covar_L
      # This is the (0,0) entry in the covariance matrix
      self.proposal_sigma = self.covar_mat[0, 0]


cdef class MISO_PairedEnd:
  """
  MISO paired-end engine.
  """
  pass
      
  
  
