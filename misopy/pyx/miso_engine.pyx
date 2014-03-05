##
## MISO engine
##
## Yarden Katz <yarden@mit.edu>
##
cimport cython

cimport matrix_utils
cimport array_utils
cimport math_utils
cimport sampling_utils

cdef class SingleEndEngine:
  """
  MISO single-end engine.
  """
  cdef char* sampler_type
  ## Gene information
  # Number of isoforms
  cdef int num_iso
  # Isoform lengths
  cdef int[:] iso_lens 
  
  ## Reads information
  cdef public int[:, :] reads
  cdef int read_len
  cdef int overhang_len

  # Sampler parameters
  cdef int num_iters
  cdef int burn_in
  cdef int num_chains
  cdef int lag
  # Number of samples to take 
  cdef int num_samples 
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
      self.sampler_type = "single-end"
      
  
  def __init__(self, reads, iso_lens, read_len, covar_mat, covar_L,
               psi_samples,
               overhang_len=1,
               num_iters=2500,
               burn_in=500,
               num_chains=1,
               lag=10):
      """
      Initialize parameters for Python instantiation of object.
      """
      self.reads = reads
      self.read_len = read_len
      self.num_iso = iso_lens.shape[0]
      self.covar_mat = covar_mat
      self.covar_L = covar_L
      # Where Psi samples will be stored (2d double array)
      self.psi_samples = psi_samples
      self.overhang_len = overhang_len
      self.num_iters = num_iters
      self.burn_in = burn_in
      self.num_chains = num_chains
      self.lag = lag
      # This is the (0,0) entry in the covariance matrix
      self.proposal_sigma = self.covar_mat[0, 0]
      # Calculate the number of samples to generate, 
      # incorporating burn-in and lag
      self.num_samples = \
        <int>((self.num_iters - self.burn_in) / self.lag) 
      # Multiply by the number of chains
      self.num_samples = self.num_chains * self.num_samples

      
  def filter_reads(self):
      """
      Filter reads.
      """
      pass

      
  def run_sampler(self):
      """
      Run sampler. This is the main Metropolis-Hastings loop:

      (1) Initialize Psi values and assignments (uniformly).

      (2) Propose Psi_next and accept with MH ratio.

      (3) For each read, sample reassignment to one of the isoforms.
      """
      cdef int n_iter = 0
      cdef int lag_counter = 0
      # Number of kept samples
      cdef int kept_samples = 0
      # Whether or not to track lag
      cdef int track_lag = 0
      # Proposals statistics
      cdef int num_rejected = 0
      cdef int num_accepted = 0
      # Initial values for psi and alpha vectors
      cdef double[:] init_psi_vector = \
        array_utils.get_double_array(self.num_iso)
      # Initialize psi and alpha vectors to be 1/k (k many isoforms)
      init_psi_vector[:] = 1/(<double>self.num_iso)
      cdef double[:] init_alpha_vector = \
        array_utils.get_double_array(self.num_iso - 1)
      init_alpha_vector[:] = 1/(<double>self.num_iso)
      # Log scores to return
      cdef double[:] log_scores = \
        array_utils.get_double_array(self.num_samples)
      print "Running %s" %(self.sampler_type)
      # First filter reads
      self.filter_reads()
      for n_iter in xrange(self.num_iters):
          if (n_iter == (self.burn_in - 1)):
              print "Counting %d as start of post burn in" %(n_iter)
              track_lag = 1
          # Calculate the MH ratio here
          # ...

          # miso_proposals.compute_mh_ratio(reads, assignments,
          #                                 new_psi_vector, new_alpha_vector,
          #                                 curr_psi_vector, curr_alpha_vector,
          #                                 gene,
          #                                 proposal_type="norm_drift",
          #                                 hyperparams=hyperparams)
          curr_psi = 0.0
          curr_alpha = 0.0
          # Keep sample and reset counter
          if (lag_counter == self.lag):
              # Store samples and their log scores
              #psi_samples[kept_samples] = curr_psi
              #log_scores[kept_samples] = curr_log_score
              kept_samples += 1
              print "Resetting lag on %d" %(n_iter)
              lag_counter = 0
          # Advance lag if we're supposed to keep track of lag,
          # i.e. if we're past burn in
          if (track_lag == 1):
              print "Tracking lag on iteration %d" %(n_iter)
              lag_counter += 1
          # Reset keeping sample
          keep_sample = 0


#     def run_sampler(self, num_iters, reads, gene, hyperparameters, params,
#                     output_file,
#                     num_chains=1,
#                     burn_in=1000,
#                     lag=2):
#         """
#         Main Metropolis-Hastings loop:

#         (1) Initialize Psi values and assignments (uniformly).

#         (2) Propose Psi_next and accept with MH ratio.

#         (3) For each read, sample reassignment to one of the available isoforms.
#         """
#         #print >> sys.stderr, "Running on %s" %(gene.label)
#         self.num_isoforms = len(gene.isoforms)
#         num_isoforms = self.num_isoforms
#         # Record gene
#         self.gene = gene
#         if self.paired_end:
#             reads = self.filter_improbable_reads(reads)

#             # Construct a matrix of considered fragment lengths by the number
#             # of isoforms (for vectorization purposes)
#             self.frag_range_matrix = transpose(tile(self.frag_range,
#                                                     [num_isoforms, 1]))

#         if len(reads) == 0:
#             print "No reads for gene: %s" %(self.gene.label)
#             print reads
#             return

#         #output_file = output_file + ".%d_iters.%d_burnin.%d_lag" %(num_iters, burn_in, lag)
#         # Don't need to put parameters in filename
#         output_file = output_file + ".miso"
#         # If output filename exists, don't run sampler
#         if os.path.isfile(os.path.normpath(output_file)):
#             print "Output filename %s exists, not running MISO." %(output_file)
#             return None
#         self.params['iters'] = num_iters
#         self.params['burn_in'] = burn_in
#         self.params['lag'] = lag
#         # Define local variables related to reads and overhang
#         self.overhang_len = self.params['overhang_len']
#         self.read_len = self.params['read_len']
#         ##
#         ## Precompute fixed parameters related to sampling
#         ##
#         # Read length
#         self.num_reads = len(reads)
#         # Number of exons per isoform
#         self.num_parts_per_isoform = \
#             np.array([len(isoform.parts) for isoform in gene.isoforms])
#         # Length for each isoform
#         self.iso_lens = \
#             np.array([isoform.len for isoform in gene.isoforms])
#         # Isoform numbers: 0,...,K-1 for K isoforms
#         self.iso_nums = np.arange(num_isoforms)
#         # Log number of reads possible per isoform
#         num_overhang_excluded = \
#             2*(self.overhang_len - 1)*(self.num_parts_per_isoform[self.iso_nums] - 1)
#         num_reads_possible = \
#             (self.iso_lens[self.iso_nums] - self.read_len + 1) - num_overhang_excluded
#         # The number of reads possible is the number of ways a read can be
#         # aligned onto the length of the current isoforms, minus the number of positions
#         # that are ruled out due to overhang constraints.
#         self.log_num_reads_possible_per_iso = \
#             np.log((self.iso_lens[self.iso_nums] - self.read_len + 1) - num_overhang_excluded)
#         self.scaled_lens_single_end = gene.iso_lens - self.read_len + 1
#         # Record cases where read length is greater than isoform length
#         # as zero (i.e. not possible)
#         self.scaled_lens_single_end[where(self.scaled_lens_single_end <= 0)] = 0
#         ##
#         ## TODO: Do same for paired-end!
#         ##
#         self.scaled_lens_paired_end = None


#         rejected_proposals = 0
#         accepted_proposals = 0
#         psi_vectors = []
# #        log_scores = {}
#         all_psi_proposals = []
#         proposal_type = "drift"
#         init_psi = np.ones(num_isoforms)/float(num_isoforms)
#         # Initialize the starting Psi vector randomly if it's the two isoform case
#         if num_isoforms == 2:
#             init_psi = dirichlet(ones(num_isoforms)/float(num_isoforms))
#         # Do not process genes with one isoform
#         elif num_isoforms == 1:
#             one_iso_msg = "Gene %s has only one isoform; skipping..." \
#                           %(self.gene.label)
#             self.miso_logger.info(one_iso_msg)
#             self.miso_logger.error(one_iso_msg)
#             return
#         curr_psi_vector = init_psi
#         # Initialize assignments of reads to isoforms in a consistent way
#         assignments = self.initialize_valid_assignments(reads)
#         # Main sampler loop
#         lag_counter = 0
#         ##
#         ## Initial alpha vector to be all ones
#         ##
#         curr_alpha_vector = log(ones(num_isoforms - 1) * (float(1)/(num_isoforms-1)))
#         ##
#         ## Change to initial condition of alpha vector!  Make it equivalent to the initial Psi value
#         ##
#         #curr_alpha_vector = log(init_psi[:-1])
#         burn_in_counter = 0
#         total_log_scores = []
#         kept_log_scores = []
#         # Compute Metropolis ratio: use drift proposal
#         proposal_score_func = self.log_score_norm_drift_proposal
#         # Set Psi vectors from proposal
#         new_psi_vector, new_alpha_vector = \
#             self.propose_psi_vector(curr_psi_vector, curr_alpha_vector)
#         curr_psi_vector = new_psi_vector
#         curr_alpha_vector = new_alpha_vector

#         already_warned = False
#         for curr_iter in xrange(num_iters):
#             # Propose a Psi value
#             new_psi_vector, new_alpha_vector = \
#                 self.propose_psi_vector(curr_psi_vector, curr_alpha_vector)
#             all_psi_proposals.append(new_psi_vector[0])
#             if curr_iter > 0:
#                 m_ratio, curr_joint_score, proposed_joint_score = \
#                  self.compute_metropolis_ratio(reads, assignments,
#                                                new_psi_vector, new_alpha_vector,
#                                                curr_psi_vector, curr_alpha_vector,
#                                                proposal_score_func, gene,
#                                                hyperparameters)
#             else:
#                 m_ratio, curr_joint_score, proposed_joint_score = \
#                     self.compute_metropolis_ratio(reads, assignments,
#                                                   new_psi_vector, new_alpha_vector,
#                                                   curr_psi_vector, curr_alpha_vector,
#                                                   proposal_score_func, gene,
#                                                   hyperparameters,
#                                                   full_metropolis=False)
#             acceptance_prob = min(1, m_ratio)
#             if rand() < acceptance_prob:
#                 jscore = proposed_joint_score
#                 # Accept sample
#                 curr_psi_vector = new_psi_vector
#                 curr_alpha_vector = new_alpha_vector
#                 accepted_proposals += 1
#             else:
#                 jscore = curr_joint_score
#                 rejected_proposals += 1

#             # Error check the log joint score
#             if isnan(jscore):
#                 self.miso_logger.error("Unable to get log joint score for gene %s" \
#                                        %(gene.label))

#             if burn_in_counter >= burn_in:
#                 # Accumulate Psi vectors
#                 if (lag_counter == lag - 1):
#                     lag_counter = 0
#                     psi_vectors.append(curr_psi_vector)
#                     kept_log_scores.append(jscore)
#                     curr_joint_score = \
#                         self.log_score_joint(reads,
#                                              assignments,
#                                              curr_psi_vector,
#                                              gene,
#                                              hyperparameters)
#                 else:
#                     lag_counter += 1
#             else:
#                 curr_joint_score = \
#                     self.log_score_joint(reads, assignments, curr_psi_vector, gene,
#                                          hyperparameters)
#             total_log_scores.append(jscore)
#             burn_in_counter += 1

#             ##
#             ## For each read, sample its reassignment to one of the isoforms
#             ##
#             # Get the log scaled Psi (psi frag) for current Psi vector
#             curr_log_psi_frag = \
#                 miso_scores.compute_log_psi_frag(curr_psi_vector,
#                                                  self.scaled_lens_single_end,
#                                                  self.num_isoforms)
#             reassignments = self.sample_reassignments(reads, curr_log_psi_frag, gene)
#             if len(reassignments) == 0:
#                 empty_reassign_msg = "Empty reassignments for reads! " + str(reads)
#                 self.miso_logger.error(empty_reassign_msg)
#                 raise Exception, empty_reassign_msg
#             curr_joint_score = \
#                 self.log_score_joint(reads,
#                                      assignments,
#                                      curr_psi_vector,
#                                      gene,
#                                      hyperparameters)
#             if curr_joint_score == -inf:
#                 self.miso_logger.error("Moved to impossible state!")
#                 self.miso_logger.error("reassignments: " + str(reassignments))
#                 self.miso_logger.error("reads: " + str(reads))
#                 raise Exception, "Moved to impossible state!"
#             assignments = reassignments

#         if accepted_proposals == 0:
#             self.miso_logger.error("0 proposals accepted!")
#             raise Exception, "0 proposals accepted!"
#         percent_acceptance = \
#             (float(accepted_proposals)/(accepted_proposals + rejected_proposals))*100
#         # Write output to file
#         print "Outputting samples to: %s..." %(output_file)
#         classes_to_counts = reads_utils.collapse_isoform_assignments(reads)
#         read_classes = []
#         read_class_counts = []
#         for read_class in classes_to_counts:
#             read_classes.append(read_class)
#             read_class_counts.append(classes_to_counts[read_class])
#         reads_data = (read_classes, read_class_counts)
#         self.output_miso_results(output_file,
#                                  gene,
#                                  reads_data,
#                                  assignments,
#                                  psi_vectors,
#                                  kept_log_scores,
#                                  num_iters,
#                                  burn_in,
#                                  lag,
#                                  percent_acceptance,
#                                  proposal_type)
              
  
    

cdef class PairedEndEngine:
  """
  MISO paired-end engine.
  """
  pass
      
  
  
