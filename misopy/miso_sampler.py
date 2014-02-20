##
## MCMC sampler for Mixture-of-Isoforms (MISO) model
##
## Yarden Katz <yarden@mit.edu>
##
## The sampler uses a Metropolis-Hastings sampling scheme, combined with
## a Gibbs sampling step.
##
#import pygsl
#from pygsl.rng import dirichlet_lnpdf
import scipy
import misopy
import numpy
import numpy as np
from numpy.random import multinomial
from numpy.random import multivariate_normal
from numpy.random import normal
from numpy import linalg as la
import misopy.reads_utils as reads_utils
from misopy.read_simulator import simulate_reads, print_reads_summary, \
                                  read_counts_to_read_list, \
                                  get_reads_summary
from scipy.misc import logsumexp
from scipy.special import gammaln
import misopy.hypothesis_test as ht
from misopy.Gene import Gene, Exon
from misopy.py2c_gene import *

# C MISO interface
import pysplicing

# Cython interface
import misopy.miso_scores as miso_scores

from scipy import *
from numpy.random.mtrand import dirichlet
from numpy import *
import cPickle as pickle
from scipy.stats import mode
import math
import time
from numpy import numarray
import os
import sys
from collections import defaultdict
import glob
import logging
import logging.handlers

import math
import operator

loggers = {}


def dirichlet_lnpdf(alpha, x):
    """
    Substitute for dirichlet_lnpdf of pygsl.
    """
    dir_log_pdf = \
        gammaln(sum(alpha)) - sum(gammaln(alpha)) + numpy.dot((alpha - 1).T, log(x).T)
    return dir_log_pdf


def get_logger(logger_name, log_outdir,
               level=logging.WARNING,
               include_stdout=True):
    """
    Return a logging object.
    """
    global loggers
    # Avoid race-conditions
    try:
        os.makedirs(log_outdir)
    except OSError:
        pass
    if loggers.get(logger_name):
        return loggers.get(logger_name)
    logger = logging.getLogger(logger_name)
    formatter = \
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                          datefmt='%m/%d/%Y %I:%M:%S %p')
    # Do not log to file
    #if log_outdir is not None:
    #    log_filename = os.path.join(log_outdir, "%s.log" %(logger_name))
    #    fh = logging.FileHandler(log_filename)
    #    fh.setLevel(level)
    #    fh.setFormatter(formatter)
    #    logger.addHandler(fh)
    logging.root.setLevel(level)
    # Optionally add handler that streams all logs
    # to stdout
    if include_stdout:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    logger.info("Created logger %s" %(logger_name))
    loggers.update({logger_name: logger})
    return logger


##
## Helper statistics/linear algebra functions
##
def logit_inv(x):
    """
    Takes a value on x \in (-inf, inf) and transforms it to a value on (0, 1).
    """
    #p = exp(x)/(1+exp(x))
    denom = append(x, 0)
    p = exp(x)/(sum(exp(denom)))
    return p

def logit(p):
    """
    Takes a value p \in (0, 1) and transforms it to (-inf, inf).
    """
    x = log(p/(1-p))
    return x

def set_diag(a, v):
    for i, elt in enumerate(a):
        a[i, i] = v
    return a

def maxi(l):
    m = max(l)
    for i, v in enumerate(l):
        if m == v:
            return i

def mini(l):
    m = min(l)
    for i, v in enumerate(l):
        if m == v:
            return i


def exp_logsumexp(a):
    return exp(a - logsumexp(a))


def vect_logsumexp(a, axis=None):
    """
    Assume a is an array.
    """
    if axis is None:
        return logsumexp(a)
    #a = asarray(a)
    shp = list(a.shape)
    shp[axis] = 1
    a_max = a.max(axis=axis)
    s = log(exp(a - a_max.reshape(shp)).sum(axis=axis))
    lse  = a_max + s
    return lse


def print_assignment_summary(assignments):
    counts = defaultdict(int)
    for a in assignments:
        counts[a] += 1
    for k, v in counts.iteritems():
        print "Total of %d in isoform %d" %(v, k)


def float_array_to_str(array_of_floats):
    """
    Convert a float numpy array to a string for printing purposes.
    """
    str_float_array = '[' + ' '.join(['%.3f' %(val) for val in array_of_floats]) + ']'
    return str_float_array


def get_paired_end_sampler_params(num_isoforms,
                                  mean_frag_len,
                                  frag_variance,
                                  read_len,
                                  overhang_len=1):
    """
    Return parameters for MISO sampler, in paired-end mode.
    """
    hyperparameters = ones(num_isoforms)
    proposal_diag = 0.05
    sigma = set_diag(zeros([num_isoforms-1, num_isoforms-1]),
                     proposal_diag)
    sampler_params = {'read_len': read_len,
                      'overhang_len': overhang_len,
                      'uniform_proposal': False,
                      'sigma_proposal': sigma,
                      'mean_frag_len': mean_frag_len,
                      'frag_variance': frag_variance}
    return sampler_params


def get_single_end_sampler_params(num_isoforms,
                                  read_len,
                                  overhang_len=1):
    """
    Return parameters for MISO sampler, in single-end mode.
    """
    hyperparameters = ones(num_isoforms)
    proposal_diag = 0.05
    sigma = set_diag(zeros([num_isoforms-1, num_isoforms-1]),
                     proposal_diag)
    sampler_params = {'read_len': read_len,
                      'overhang_len': overhang_len,
                      'uniform_proposal': False,
                      'sigma_proposal': sigma}
    return sampler_params


class MISOSampler:
    def __init__(self, params,
                 paired_end=False,
                 log_dir=None):
        """
        Make a sampler with the given parameters.
        """
        self.params = params
        self.paired_end = paired_end
        # set default fragment length distribution parameters
        if self.paired_end:
            if ((not 'mean_frag_len' in self.params) or \
                (not 'frag_variance' in self.params)):
                raise Exception, "Must set mean_frag_len and frag_variance when " \
                      "running in sampler on paired-end data."
            self.mean_frag_len = self.params['mean_frag_len']
            self.frag_variance = self.params['frag_variance']

        if log_dir != None:
            self.log_dir = os.path.abspath(os.path.expanduser(log_dir))
            self.log_dir = os.path.join(log_dir, 'logs')
            # Avoid race-conditions
            try:
                os.makedirs(self.log_dir)
            except OSError:
                pass
        self.miso_logger = get_logger('miso_logger', self.log_dir)
        self.miso_logger.info("Instantiated sampler.")
        self.gene = None
        self.iso_lens_matrix = None


    # def log_multivar_normal_pdf(x, mu, sigma):
    #     """
    #     Multivariate Normal pdf (adapted from PyMix.)
    #     """
    #     dd = la.det(sigma)
    #     inverse = la.inv(sigma)
    #     N = len(x)
    #     ff = math.pow(2*math.pi,-N/2.0)*math.pow(dd,-0.5)
    #     # centered input values
    #     centered = numarray.subtract(x, numarray.repeat([mu], N))
    #     res = ff * numarray.exp(-0.5*numarray.sum(numarray.multiply(centered,numarray.dot(centered,inverse)), 1))
    #     return numarray.log(res)[0]


    def logistic_normal_log_pdf(self, theta, mu):
        """
        The log of the PDF for the multivariate Logistic-Normal distribution.
        Written in terms of k-1 dimensions.
        """
        theta = theta[:-1]
        if len(theta) != len(mu):
            raise Exception, "len(theta) = %d != len(mu) = %d" \
                  %(len(theta), len(mu))
        theta_last = 1-sum(theta)
        sigma = self.params['sigma_proposal']
        covar_constant = power(float(la.det(2*math.pi*sigma)), -0.5)
        invsigma = la.inv(transpose(sigma))
        prod_theta = 1/(float(prod(theta))*theta_last)
        first_log = -0.5*transpose(log(theta/theta_last) - mu)
        second_log = log(theta/theta_last) - mu
        exp_part = dot(dot(first_log, invsigma), second_log)
        pdf_value = covar_constant*prod_theta*exp(exp_part)
        return log(pdf_value)


    def log_score_joint(self, reads, assignments, psi_vector,
                        gene, hyperparameters):
        """
        Return a log score for the joint distribution.
        """
        # Get the logged Psi frag vector
        # Score the read
        if not self.paired_end:
#            log_reads_prob = \
#                sum(self.log_score_reads(reads, assignments, gene))
#           log_reads_prob = \
#                sum(miso_scores.log_score_reads(reads,
#                                                assignments,
#                                                self.num_parts_per_isoform,
#                                                self.iso_lens,
#                                                self.num_reads,
#                                                self.log_num_reads_possible_per_iso))
            log_reads_prob = \
                miso_scores.sum_log_score_reads(reads,
                                                assignments,
                                                self.num_parts_per_isoform,
                                                self.iso_lens,
                                                self.num_reads,
                                                self.log_num_reads_possible_per_iso)
            log_psi_frag = \
                miso_scores.compute_log_psi_frag(psi_vector,
                                                 self.scaled_lens_single_end,
                                                 self.num_isoforms)
        else:
            raise Exception, "Paired-end not implemented!"
            log_reads_prob = \
                sum(self.log_score_paired_end_reads(reads, assignments, gene))
            log_psi_frag = None
        if not self.paired_end:
            log_assignments_prob = \
                miso_scores.sum_log_score_assignments(assignments,
                                                      log_psi_frag,
                                                      self.num_reads)
        else:
            raise Exception, "Paired-end not implemented!"
            log_assignments_prob = \
                sum(self.log_score_paired_end_assignment(assignments,
                                                         psi_vector,
                                                         gene))
        # Score the Psi vector: keep this in Python for now
        #log_psi_prob = self.log_score_psi_vector(psi_vector, hyperparameters)
        log_psi_prob = \
            miso_scores.log_score_psi_vector(psi_vector, hyperparameters)
        log_joint_score = log_reads_prob + log_assignments_prob + log_psi_prob
        return log_joint_score

#     def log_score_joint(self, reads, assignments, psi_vector, gene, hyperparameters):
#         """
#         Return a log score for the joint distribution.
#         """
#         log_reads_prob = 0
#         log_assignments_prob = 0
#         #for read, read_assignment in zip(reads, assignments):
#         for n in xrange(len(reads)):
#             read = reads[n]
#             read_assignment = assignments[n]
#             # Score all the reads
#             #log_reads_prob += self.log_score_read(read, read_assignment, gene)
#             # Score reads using table
#             log_reads_prob += self.log_reads_table[(read, read_assignment)]
#             # Score all the assignments
#             log_assignments_prob += self.log_score_assignment(read_assignment, psi_vector,
#                                                               gene)
#         # Score the Psi vector
#         log_psi_prob = self.log_score_psi_vector(psi_vector, hyperparameters)
#         log_joint_score = log_reads_prob + log_assignments_prob + log_psi_prob
#         return log_joint_score

    def log_score_psi_vector(self, psi_vector, hyperparameters):
        """
        Score a setting of the Psi values given hyperparameters of
        the Dirichlet prior.
        """
        #assert (all(hyperparameters > 0))
        # Hyperparameters first, Psi vector second (unlike C interface)
        return dirichlet_lnpdf(hyperparameters, psi_vector)


#    def log_score_psi_vector_proposal(self, psi_vector, alpha_vector):
#        """
#        Just like log_score_psi_vector, but for a proposal distribution.
#        To keep consistent with the drift scoring function, this function
#        takes an alpha_vector argument but does not use it in the computation.
#        """
#        hyperparameters = [1]*len(psi_vector)
#        return dirichlet_lnpdf(hyperparameters, psi_vector)


    def log_score_paired_end_assignment(self, isoform_nums, psi_vector, gene):
        """
        Score an assignment of a set of paired-end reads given psi
        and a gene (i.e. a set of isoforms).

        Score each assignment by the insert length and its probability.
        """
        ###
        ### TODO: Compute this once before starting sampling!!!
        ###
        # Construct matrix of number of considered fragment lengths
        # by isoform lengths (for vectorization purposes)
        iso_lens_matrix = tile(gene.iso_lens,
                               [self.frag_range_len, 1])

        # Scaled lengths of each isoform by the considered
        # fragment length in the distribution
        scaled_lens = iso_lens_matrix - self.frag_range_matrix + 1

        # Construct a matrix of Psi values with dimensions of
        # number of considered fragment lengths by number of
        # isoforms (for vectorization purposes)
        #self.psi_matrix = tile(psi_vector,
        #                       [self.frag_range_len, 1])

        # Take the log of the valid Psi frags
        psi_frag = log(psi_vector) + log(scaled_lens)

        # Convert NaN to -Inf
        psi_frag[where(isnan(psi_frag))[0]] = -Inf

        # Create a masked array where the elements that are -Inf
        masked_psi_frag = ma.masked_where(psi_frag == -Inf,
                                          psi_frag)

        # Normalize the scaled Psi for each possible fragment length
        # (based on the fragment length distribution)
        psi_frag = vect_logsumexp(psi_frag, axis=0)

        psi_frag = psi_frag - logsumexp(psi_frag)

        psi_frags = tile(psi_frag, [self.num_reads, 1])
        final_psi_frags = psi_frags[range(self.num_reads),
                                    isoform_nums]
        return final_psi_frags


    # def log_score_reads(self, reads, isoform_nums, gene):
    #     """
    #     Score a set of reads given their isoform assignments.
    #     Vectorized version.
    #     """
    #     # The probability of a read being assigned to an isoform that
    #     # could not have generated it (i.e. where the read is not a
    #     # substring of the isoform) is zero.  Check for consistency
    #     overhang_excluded = \
    #         2*(self.overhang_len - 1)*(gene.num_parts_per_isoform[isoform_nums] - 1)
    #     # The number of reads possible is the number of ways a 36 nt long read can be
    #     # aligned onto the length of the current isoforms, minus the number of positions
    #     # that are ruled out due to overhang constraints.
    #     num_reads_possible = \
    #         (gene.iso_lens[isoform_nums] - self.read_len + 1) - overhang_excluded
    #     log_prob_reads = log(1) - log(num_reads_possible)
    #     zero_prob_indx = nonzero(reads[np.arange(self.num_reads), isoform_nums] == 0)[0]
    #     # Assign probability 0 to reads inconsistent with assignment
    #     log_prob_reads[zero_prob_indx] = -inf
    #     return log_prob_reads


    def log_score_paired_end_reads(self, reads, isoform_nums, gene):
        """
        Score a set of paired-end reads given their isoform assignments.
        Vectorized version.

        Takes in a set of paired-end reads (set of pairs of read alignments and
        their fragmentlengths), a set of isoform assignments for those reads
        (isoforms by their number), and a gene.

        Optional arguments:

          - frag_len_dist: a function that takes in a fragment length and
            returns a log probability. By default, this is set to a binomial
            distribution (with a mean-variance parameterization)
        """
        # The probability of a read being assigned to an isoform that
        # could not have generated it (i.e. where the read is not a
        # substring of the isoform) is zero.
        pe_reads = reads[:, 0]
        frag_lens = reads[:, 1]

        # Get the number of overhang positions violated
        num_overhang_excluded = 0
        # The number of reads possible is the number of ways a fragment of the given length
        # can be generated from the isoforms
#       print "isoform_nums: ", isoform_nums
        assigned_iso_frag_lens = frag_lens[range(self.num_reads), isoform_nums]
#        print "assigned_iso_frag_lens: ", assigned_iso_frag_lens
#       print >> sys.stderr, "Scoring fragments: ", assigned_iso_frag_lens, " - reads: ", pe_reads
#       print >> sys.stderr, "Isoform assignments: ", isoform_nums
        num_reads_possible = \
            gene.iso_lens[isoform_nums] - assigned_iso_frag_lens + 1 - num_overhang_excluded
#       print >> sys.stderr, "num_reads_possible: ", num_reads_possible
        # The probability of a paired-end read is
        # (1 / num_possible_reads) * P(fragment_len)
#       print >> sys.stderr, "log_frag_len_prob: ", self.log_frag_len_prob(assigned_iso_frag_lens)
#       print >> sys.stderr, " reads : ", type(pe_reads)
        log_prob_frags = self.log_frag_len_prob(assigned_iso_frag_lens,
                                                self.mean_frag_len, self.frag_variance)
        log_prob_reads = (log(1) - log(num_reads_possible)) + log_prob_frags
#        print "log_prob_reads: ", log_prob_reads
        zero_prob_indx = nonzero(pe_reads[range(self.num_reads), isoform_nums] == 0)[0]
        # Assign probability 0 to reads inconsistent with assignment
        log_prob_reads[zero_prob_indx] = -inf
        return log_prob_reads


    # def log_score_assignment(self, isoform_nums, psi_vector):
    #     """
    #     Score an assignment of a set of reads given psi
    #     and a gene (i.e. a set of isoforms).
    #     """
    #     ###
    #     ### TODO: no need to compute scaled lens each time!
    #     ###
    #     psi_frag = log(psi_vector) + log(self.scaled_lens_single_end)
    #     psi_frag = psi_frag - logsumexp(psi_frag)
    #     psi_frags = tile(psi_frag, [self.num_reads, 1])
    #     # NEW VERSION: uses xrange
    #     return psi_frags[np.arange(self.num_reads), isoform_nums]
    #     # OLD VERSION: uses range
    #     #return psi_frags[range(self.num_reads), isoform_nums]


    def propose_psi_vector(self, psi_vector, alpha_vector):
        """
        Propose a new Psi vector.
        """
        # Independent uniform proposal distribution
        proposed_psi_vector = None
        proposed_alpha_vector = []
        #if len(alpha_vector) != 0:
        proposed_psi_vector, proposed_alpha_vector = \
            self.propose_norm_drift_psi_alpha(alpha_vector)
        proposed_psi_vector = np.array(proposed_psi_vector)
        return (proposed_psi_vector, proposed_alpha_vector)


    def propose_norm_drift_psi_alpha(self, alpha_vector):
        if len(alpha_vector) == 1:
            alpha_vector = alpha_vector[0]
#            print "proposing from normal with mean: ", alpha_vector, " exp: ", exp(alpha_vector)
            alpha_next = [normal(alpha_vector, self.params['sigma_proposal'])]
#            print "got alpha_next: ", alpha_next, " exp: ", exp(alpha_next)
            new_psi = logit_inv([alpha_next[0]])[0]
            new_psi_vector = [new_psi, 1-new_psi]
        else:
            alpha_next = multivariate_normal(alpha_vector, self.params['sigma_proposal'])
            new_psi = logit_inv(alpha_next)
            new_psi_vector = concatenate((new_psi, array([1-sum(new_psi)])))
        return (new_psi_vector, alpha_next)


    def log_score_norm_drift_proposal(self, target_psi, alpha_vector,
                                      hyperparameters=[]):
        """
        Compute the log probability of transitioning to the target_psi given
        that the proposal distribution has a mean of source_alpha_vector.
        """
        log_transition = \
            self.logistic_normal_log_pdf(target_psi, alpha_vector)
        return log_transition


    def sample_reassignments(self, reads, log_psi_frag, gene):
        """
        Sample a reassignments of reads to isoforms.
        Note that this does not dependent on the read's current assignment since
        we're already considering the possibility of 'reassigning' the read to
        its current assignment in the probability calculations.
        """
        reassignment_probs = []
        all_assignments = transpose(tile(arange(self.num_isoforms, dtype=int),
                                         [self.num_reads, 1]))
        for assignment in all_assignments:
            if not self.paired_end:
                # Single-end
                # Score reads given their assignment
                read_probs = \
                    miso_scores.log_score_reads(reads,
                                                assignment,
                                                self.num_parts_per_isoform,
                                                self.iso_lens,
                                                self.num_reads,
                                                self.log_num_reads_possible_per_iso)
                # Score the assignments of reads given Psi vector
                assignment_probs = \
                    miso_scores.log_score_assignments(assignment,
                                                      log_psi_frag,
                                                      self.num_reads)
            else:
                raise Exception, "Paired-end unimplemented!"
                # Paired-end
                # Score paired-end reads given their assignment
                read_probs = \
                    self.log_score_paired_end_reads(reads, assignment, gene)
                # Score the assignments of paired-end reads given Psi vector
                #assignment_probs = \
                #    self.log_score_paired_end_assignment(assignment, psi_vector, gene)
            reassignment_p = read_probs + assignment_probs
            reassignment_probs.append(reassignment_p)
        reassignment_probs = transpose(array(reassignment_probs))
        m = transpose(vect_logsumexp(reassignment_probs, axis=1)[newaxis,:])
        norm_reassignment_probs = exp(reassignment_probs - m)

        rvsunif = random.rand(self.num_reads, 1)
        yrvs = (rvsunif<cumsum(norm_reassignment_probs,axis=1)).argmax(1)[:,newaxis]
        ### Note taking first element of transpose(yrvs)!  To avoid a list of assignments
        return transpose(yrvs)[0]


    def compute_metropolis_ratio(self, reads, assignments, proposed_psi_vector,
                                 proposed_alpha_vector,
                                 curr_psi_vector, curr_alpha_vector,
                                 proposal_score_func, gene,
                                 hyperparameters,
                                 full_metropolis=True):
        """
        Compute the Metropolis-Hastings ratio:

            P(psi_next)Q(psi; psi_next)
            ---------------------------
            P(psi)Q(psi_next; psi)
        """
        # Compute acceptance ratio: the joint score for proposed Psi divided
        # by joint score given current Psi
        # P(Psi', ...)
        proposed_joint_score = \
            self.log_score_joint(reads, assignments, proposed_psi_vector,
                                 gene, hyperparameters)
        # P(Psi, ...)
        curr_joint_score = \
            self.log_score_joint(reads, assignments, curr_psi_vector,
                                 gene, hyperparameters)
        if curr_joint_score == -inf:
            self.miso_logger.error("Joint score of current state is negative infinity!")
            self.miso_logger.error("  - assignments: " + str(assignments))
            self.miso_logger.error("  - psi vector: " + str(curr_psi_vector))
            self.miso_logger.error("  - reads: " + str(reads))
            raise Exception, "curr_joint_score is negative."
        # Q(x; x'), the probability of proposing to move back to current state from
        # proposed state x'
        mh_ratio = None
        proposal_to_curr_score = \
            proposal_score_func(curr_psi_vector, proposed_alpha_vector)
        # Q(x'; x), the probability of proposing to move to the proposed state x' from
        # the current state
        curr_to_proposal_score = \
            proposal_score_func(proposed_psi_vector, curr_alpha_vector)
#       self.miso_logger.debug("Computing MH ratio...")
#       self.miso_logger.debug("  - Proposed Psi vector: " + str(proposed_psi_vector))
#       self.miso_logger.debug("  - curr_joint_score: " + str(curr_joint_score))
#       self.miso_logger.debug("  - proposed_joint_score: " + str(proposed_joint_score))
#       self.miso_logger.debug("  - curr_to_proposal_score: " + str(curr_to_proposal_score))
#       self.miso_logger.debug("  - proposal_to_curr_score: " + str(proposal_to_curr_score))
#       self.miso_logger.debug("  - assignments: " + str(assignments))

#       if abs(curr_to_proposal_score - proposal_to_curr_score) > 2:
#           self.miso_logger.warn("curr_to_proposal and proposal_to_curr diverge!")
#           self.miso_logger.warn("proposed_psi_vector: " + str(proposed_psi_vector))
#           self.miso_logger.warn("proposed_alpha_vector: " + str(proposed_alpha_vector))
#           self.miso_logger.warn("curr_psi_vector: " + str(curr_psi_vector))
#           self.miso_logger.warn("curr_alpha_vector: " + str(curr_alpha_vector))
#           time.sleep(4)
#           raise Exception
        # Computing full Metropolis-Hastings ratio
        if not full_metropolis:
            mh_ratio = (proposed_joint_score - curr_joint_score)
        else:
            mh_ratio = (proposed_joint_score + proposal_to_curr_score) - \
                       (curr_joint_score + curr_to_proposal_score)
        if curr_to_proposal_score == -inf:
            self.miso_logger.error("curr to proposal is -inf")
            raise Exception, "curr to proposal is -Inf"
        if proposed_joint_score == -inf:
            self.miso_logger.debug("Proposing to move to impossible state!")
            raise Exception, "Proposing to move to impossible state!"
        if abs(mh_ratio) == Inf:
            self.miso_logger.debug("MH ratio is Inf!")
            raise Exception, "MH ratio is Inf!"
        return (exp(mh_ratio), curr_joint_score, proposed_joint_score)


    def choose_assignment(self, prob_assignments, valid_assignments):
        """
        Choose an assignment from a normalized vector of probabilities and
        a set of valid assignments.
        """
#       rand_indx = list(multinomial(1, prob_assignments)).index(1)

        assignments = list(multinomial(1, prob_assignments))

        # Ensure against bad assignments
        if 1 not in assignments:
            self.miso_logger.error("choose_assignment cannot find proper assignment.")
            raise Exception, "choose_assignment: Proper assignment cannot be found."

        # Choose assignment
        rand_indx = assignments.index(1)

        valid_assignment_indx = valid_assignments[rand_indx]
        return valid_assignment_indx


    def initialize_valid_assignments(self, reads):
        """
        Pick a random but valid/consistent assignments for the read to start with.
        """
        frag_lens = []
        if self.paired_end:
            # if paired-end, consider the alignments of the reads
            pe_reads = reads[:, 0]
            frag_lens = reads[:, 1]
        assignments = []
        if self.paired_end:
            for r, frags in zip(pe_reads, frag_lens):
                valid_assignments = nonzero(array(r))[0]
                # for paired end reads, weigh the valid assignments by the prior probability
                # of the fragment lengths they posit
                prob_assignments = []
                assignment_frag_lens = frags[valid_assignments]
                if self.paired_end:
                    for v, frag_len in zip(valid_assignments, assignment_frag_lens):
                        # Compute score of fragment length
                        frag_len_score = exp(self.log_frag_len_prob(frag_len,
                                                                    self.mean_frag_len,
                                                                    self.frag_variance))

                        # Initially, avoid assignments where reads are assigned
                        # to isoforms
                        if frag_len < self.mean_frag_len:
                            frag_len_score = 0
                        prob_assignments.append(frag_len_score)
                    prob_assignments = array(prob_assignments)
                    # renormalize
                    prob_assignments = prob_assignments / sum(prob_assignments)
                if all(isnan(prob_assignments)):
                    raise Exception, "Cannot find valid assignment for read"
                chosen_assignment = self.choose_assignment(prob_assignments, valid_assignments)
                assignments.append(chosen_assignment)
        else:
            for r in reads:
                valid_assignments = nonzero(array(r))[0]
                ###
                ### TODO: Choose assignments probabilistically WITHOUT multinomial! Use randint
                ###
                #prob_assignments = [1/float(len(valid_assignments)) for v in valid_assignments]
                #chosen_assignment = self.choose_assignment(prob_assignments, valid_assignments)
                # Choose first assignment
                chosen_assignment = valid_assignments[0]
                assignments.append(chosen_assignment)
        assignments = array(assignments, dtype=int)
        return assignments


    def filter_improbable_reads(self, reads):
        """
        Filter improbable reads based on fragment lengths, i.e.
        reads that are probabilistically inconsistent with all
        isoforms.
        """
        if len(reads) == 0:
            return array([])
        pe_reads = reads[:, 0]
        frag_lens = reads[:, 1]
        filtered_reads = []
        num_skipped = 0

        for r, frags in zip(pe_reads, frag_lens):
            valid_assignments = nonzero(array(r))[0]
            assignment_frag_lens = frags[valid_assignments]
            len_scores = []

            for v, frag_len in zip(valid_assignments,
                                   assignment_frag_lens):
                # compute the prior probability of the fragment lengths each
                # read posits for each isoform
                frag_len_score = \
                    exp(self.log_frag_len_prob(frag_len, self.mean_frag_len,
                                               self.frag_variance))
                len_scores.append(frag_len_score)

            len_scores = array(len_scores)

            # If the read posits improbable fragment lengths for each isoform,
            # discard it
            if all(len_scores == 0):
                num_skipped += 1
                continue
            else:
                filtered_reads.append([r, frags])

        print "Filtered out %d reads that posited improbable fragment lengths " \
              "with all isoforms." %(num_skipped)

        filtered_reads = array(filtered_reads)
        return filtered_reads


    # def run_sampler(self, num_iters, reads, gene, hyperparameters, params,
    #                 output_file,
    #                 num_chains=6,
    #                 burn_in=1000,
    #                 lag=2,
    #                 prior_params=None,
    #                 # By default, use sampler with read classes (collapsed)
    #                 # to get speed boost for single-end reads
    #                 # (To revert to old reassigning sampler, use
    #                 # pysplicing.MISO_ALGO_REASSIGN)
    #                 algorithm=pysplicing.MISO_ALGO_CLASSES,
    #                 start_cond=pysplicing.MISO_START_AUTO,
    #                 stop_cond=pysplicing.MISO_STOP_FIXEDNO,
    #                 verbose=True):
    #     """
    #     Fast version of MISO MCMC sampler.

    #     Calls C version and returns results.
    #     """
    #     num_isoforms = len(gene.isoforms)
    #     self.num_isoforms = num_isoforms

    #     if prior_params == None:
    #         prior_params = (1.0,) * num_isoforms

    #     read_positions = reads[0]
    #     read_cigars = reads[1]

    #     self.num_reads = len(read_positions)

    #     if self.num_reads == 0:
    #         print "No reads for gene: %s" %(gene.label)
    #         return

    #     output_file = output_file + ".miso"
    #     # If output filename exists, don't run sampler
    #     if os.path.isfile(os.path.normpath(output_file)):
    #         print "Output filename %s exists, not running MISO." \
    #             %(output_file)
    #         return None
    #     self.params['iters'] = num_iters
    #     self.params['burn_in'] = burn_in
    #     self.params['lag'] = lag

    #     # Define local variables related to reads and overhang
    #     self.overhang_len = self.params['overhang_len']
    #     self.read_len = self.params['read_len']

    #     t1 = 0
    #     t2 = 0
    #     if verbose:
    #         t1 = time.time()
    #     #self.miso_logger.info("Running sampler...")
    #     #self.miso_logger.info("  - num_iters: " + str(num_iters))
    #     #self.miso_logger.info("  - burn-in: " + str(burn_in))
    #     #self.miso_logger.info("  - lag: " + str(lag))
    #     #self.miso_logger.info("  - paired-end? " + str(self.paired_end))
    #     #self.miso_logger.info("  - gene: " + str(gene))
    #     rejected_proposals = 0
    #     accepted_proposals = 0
    #     psi_vectors = []
    #     all_psi_proposals = []

    #     if params['uniform_proposal']:
    #         self.miso_logger.debug("UNIFORM independent proposal being used.")
    #         proposal_type = "unif"
    #     else:
    #         self.miso_logger.debug("Non-uniform proposal being used.")
    #         self.miso_logger.debug("  - sigma_proposal: " + str(params['sigma_proposal']))
    #         proposal_type = "drift"
    #     init_psi = ones(num_isoforms)/float(num_isoforms)
    #     # Do not process genes with one isoform
    #     if num_isoforms == 1:
    #         one_iso_msg = "Gene %s has only one isoform; skipping..." \
    #                       %(gene.label)
    #         self.miso_logger.warning(one_iso_msg)
    #         return

    #     # Convert Python Gene object to C
    #     c_gene = py2c_gene(gene)

    #     ##
    #     ## Run C MISO
    #     ##
    #     read_positions = tuple([r+1 for r in read_positions])
    #     if self.paired_end:
    #         # Number of standard deviations in insert length
    #         # distribution to consider when assigning reads
    #         # to isoforms
    #         num_sds = 4L

    #         # Run paired-end
    #         miso_results = pysplicing.MISOPaired(c_gene, 0L,
    #                                              read_positions,
    #                                              read_cigars,
    #                                              long(self.read_len),
    #                                              float(self.mean_frag_len),
    #                                              float(self.frag_variance),
    #                                              float(num_sds),
    #                                              long(num_iters),
    #                                              long(burn_in),
    #                                              long(lag),
    #                                              prior_params,
    #                                              long(self.overhang_len),
    #                                              long(num_chains),
    #                                              start_cond,
    #                                              stop_cond)
    #     else:
    #         print "Running old sampler..."

    #         # Run single-end
    #         # miso_results = pysplicing.MISO(c_gene,
    #         #                                0L,
    #         #                                read_positions,
    #         #                                read_cigars,
    #         #                                long(self.read_len),
    #         #                                long(num_iters),
    #         #                                long(burn_in),
    #         #                                long(lag),
    #         #                                prior_params,
    #         #                                long(self.overhang_len),
    #         #                                long(num_chains),
    #         #                                start_cond,
    #         #                                stop_cond,
    #         #                                algorithm)

    #     # Psi samples
    #     psi_vectors = transpose(array(miso_results[0]))

    #     # Log scores of accepted samples
    #     kept_log_scores = transpose(array(miso_results[1]))

    #     # Read classes
    #     read_classes = miso_results[2]

    #     # Read class statistics
    #     read_class_data = miso_results[3]

    #     # Assignments of reads to isoforms
    #     assignments = miso_results[4]

    #     # Statistics and parameters about sampler run
    #     run_stats = miso_results[5]

    #     # Assignments of reads to classes.
    #     # read_classes[n] represents the read class that has
    #     # read_assignments[n]-many reads.
    #     reads_data = (read_classes, read_class_data)

    #     assignments = array(assignments)

    #     # Skip events where all reads are incompatible with the annotation;
    #     # do not output a file for those.
    #     if all(assignments == -1):
    #         print "All reads incompatible with annotation, skipping..."
    #         return

    #     accepted_proposals = run_stats[4]
    #     rejected_proposals = run_stats[5]

    #     percent_acceptance = (float(accepted_proposals)/(accepted_proposals + \
    #                                                      rejected_proposals)) * 100
    #     #self.miso_logger.info("Percent acceptance (including burn-in): %.4f" %(percent_acceptance))
    #     #self.miso_logger.info("Number of iterations recorded: %d" %(len(psi_vectors)))

    #     # Write MISO output to file
    #     print "Outputting samples to: %s..." %(output_file)
    #     self.miso_logger.info("Outputting samples to: %s" %(output_file))
    #     self.output_miso_results(output_file, gene, reads_data, assignments,
    #                              psi_vectors, kept_log_scores, num_iters,
    #                              burn_in, lag, percent_acceptance,
    #                              proposal_type)
    #     if verbose:
    #         t2 = time.time()
    #         print "Event took %.2f seconds" %(t2 - t1)



    #                 output_file,
    #                 num_chains=6,
    #                 burn_in=1000,
    #                 lag=2,

    def run_sampler(self, num_iters, reads, gene, hyperparameters, params,
                    output_file,
                    num_chains=1,
                    burn_in=1000,
                    lag=2):
        """
        Main Metropolis-Hastings loop:

        (1) Initialize Psi values and assignments (uniformly).

        (2) Propose Psi_next and accept with MH ratio.

        (3) For each read, sample reassignment to one of the available isoforms.
        """
        #print >> sys.stderr, "Running on %s" %(gene.label)
        self.num_isoforms = len(gene.isoforms)
        num_isoforms = self.num_isoforms
        # Record gene
        self.gene = gene
        if self.paired_end:
            reads = self.filter_improbable_reads(reads)

            # Construct a matrix of considered fragment lengths by the number
            # of isoforms (for vectorization purposes)
            self.frag_range_matrix = transpose(tile(self.frag_range,
                                                    [num_isoforms, 1]))

        if len(reads) == 0:
            print "No reads for gene: %s" %(self.gene.label)
            print reads
            return

        #output_file = output_file + ".%d_iters.%d_burnin.%d_lag" %(num_iters, burn_in, lag)
        # Don't need to put parameters in filename
        output_file = output_file + ".miso"
        # If output filename exists, don't run sampler
        if os.path.isfile(os.path.normpath(output_file)):
            print "Output filename %s exists, not running MISO." %(output_file)
            return None
        self.params['iters'] = num_iters
        self.params['burn_in'] = burn_in
        self.params['lag'] = lag
        # Define local variables related to reads and overhang
        self.overhang_len = self.params['overhang_len']
        self.read_len = self.params['read_len']
        ##
        ## Precompute fixed parameters related to sampling
        ##
        # Read length
        self.num_reads = len(reads)
        # Number of exons per isoform
        self.num_parts_per_isoform = \
            np.array([len(isoform.parts) for isoform in gene.isoforms])
        # Length for each isoform
        self.iso_lens = \
            np.array([isoform.len for isoform in gene.isoforms])
        # Isoform numbers: 0,...,K-1 for K isoforms
        self.iso_nums = np.arange(num_isoforms)
        # Log number of reads possible per isoform
        num_overhang_excluded = \
            2*(self.overhang_len - 1)*(self.num_parts_per_isoform[self.iso_nums] - 1)
        num_reads_possible = \
            (self.iso_lens[self.iso_nums] - self.read_len + 1) - num_overhang_excluded
        # The number of reads possible is the number of ways a read can be
        # aligned onto the length of the current isoforms, minus the number of positions
        # that are ruled out due to overhang constraints.
        self.log_num_reads_possible_per_iso = \
            np.log((self.iso_lens[self.iso_nums] - self.read_len + 1) - num_overhang_excluded)
        self.scaled_lens_single_end = gene.iso_lens - self.read_len + 1
        # Record cases where read length is greater than isoform length
        # as zero (i.e. not possible)
        self.scaled_lens_single_end[where(self.scaled_lens_single_end <= 0)] = 0
        ##
        ## TODO: Do same for paired-end!
        ##
        self.scaled_lens_paired_end = None


        rejected_proposals = 0
        accepted_proposals = 0
        psi_vectors = []
#        log_scores = {}
        all_psi_proposals = []
        proposal_type = "drift"
        init_psi = np.ones(num_isoforms)/float(num_isoforms)
        # Initialize the starting Psi vector randomly if it's the two isoform case
        if num_isoforms == 2:
            init_psi = dirichlet(ones(num_isoforms)/float(num_isoforms))
        # Do not process genes with one isoform
        elif num_isoforms == 1:
            one_iso_msg = "Gene %s has only one isoform; skipping..." \
                          %(self.gene.label)
            self.miso_logger.info(one_iso_msg)
            self.miso_logger.error(one_iso_msg)
            return
        curr_psi_vector = init_psi
        # Initialize assignments of reads to isoforms in a consistent way
        assignments = self.initialize_valid_assignments(reads)
        # Main sampler loop
        lag_counter = 0
        ##
        ## Initial alpha vector to be all ones
        ##
        curr_alpha_vector = log(ones(num_isoforms - 1) * (float(1)/(num_isoforms-1)))
        ##
        ## Change to initial condition of alpha vector!  Make it equivalent to the initial Psi value
        ##
        #curr_alpha_vector = log(init_psi[:-1])
        burn_in_counter = 0
        total_log_scores = []
        kept_log_scores = []
        # Compute Metropolis ratio: use drift proposal
        proposal_score_func = self.log_score_norm_drift_proposal
        # Set Psi vectors from proposal
        new_psi_vector, new_alpha_vector = \
            self.propose_psi_vector(curr_psi_vector, curr_alpha_vector)
        curr_psi_vector = new_psi_vector
        curr_alpha_vector = new_alpha_vector

        already_warned = False
        for curr_iter in xrange(num_iters):
            # Propose a Psi value
            new_psi_vector, new_alpha_vector = \
                self.propose_psi_vector(curr_psi_vector, curr_alpha_vector)
            all_psi_proposals.append(new_psi_vector[0])
            if curr_iter > 0:
                m_ratio, curr_joint_score, proposed_joint_score = \
                 self.compute_metropolis_ratio(reads, assignments,
                                               new_psi_vector, new_alpha_vector,
                                               curr_psi_vector, curr_alpha_vector,
                                               proposal_score_func, gene,
                                               hyperparameters)
            else:
                m_ratio, curr_joint_score, proposed_joint_score = \
                    self.compute_metropolis_ratio(reads, assignments,
                                                  new_psi_vector, new_alpha_vector,
                                                  curr_psi_vector, curr_alpha_vector,
                                                  proposal_score_func, gene,
                                                  hyperparameters,
                                                  full_metropolis=False)
            acceptance_prob = min(1, m_ratio)
            if rand() < acceptance_prob:
                jscore = proposed_joint_score
                # Accept sample
                curr_psi_vector = new_psi_vector
                curr_alpha_vector = new_alpha_vector
                accepted_proposals += 1
            else:
                jscore = curr_joint_score
                rejected_proposals += 1

            # Error check the log joint score
            if isnan(jscore):
                self.miso_logger.error("Unable to get log joint score for gene %s" \
                                       %(gene.label))

            if burn_in_counter >= burn_in:
                # Accumulate Psi vectors
                if (lag_counter == lag - 1):
                    lag_counter = 0
                    psi_vectors.append(curr_psi_vector)
                    kept_log_scores.append(jscore)
                    curr_joint_score = \
                        self.log_score_joint(reads,
                                             assignments,
                                             curr_psi_vector,
                                             gene,
                                             hyperparameters)
                else:
                    lag_counter += 1
            else:
                curr_joint_score = \
                    self.log_score_joint(reads, assignments, curr_psi_vector, gene,
                                         hyperparameters)
            total_log_scores.append(jscore)
            burn_in_counter += 1

            ##
            ## For each read, sample its reassignment to one of the isoforms
            ##
            # Get the log scaled Psi (psi frag) for current Psi vector
            curr_log_psi_frag = \
                miso_scores.compute_log_psi_frag(curr_psi_vector,
                                                 self.scaled_lens_single_end,
                                                 self.num_isoforms)
            reassignments = self.sample_reassignments(reads, curr_log_psi_frag, gene)
            if len(reassignments) == 0:
                empty_reassign_msg = "Empty reassignments for reads! " + str(reads)
                self.miso_logger.error(empty_reassign_msg)
                raise Exception, empty_reassign_msg
            curr_joint_score = \
                self.log_score_joint(reads,
                                     assignments,
                                     curr_psi_vector,
                                     gene,
                                     hyperparameters)
            if curr_joint_score == -inf:
                self.miso_logger.error("Moved to impossible state!")
                self.miso_logger.error("reassignments: " + str(reassignments))
                self.miso_logger.error("reads: " + str(reads))
                raise Exception, "Moved to impossible state!"
            assignments = reassignments

        if accepted_proposals == 0:
            self.miso_logger.error("0 proposals accepted!")
            raise Exception, "0 proposals accepted!"
        percent_acceptance = \
            (float(accepted_proposals)/(accepted_proposals + rejected_proposals))*100
        # Write output to file
        print "Outputting samples to: %s..." %(output_file)
        classes_to_counts = reads_utils.collapse_isoform_assignments(reads)
        read_classes = []
        read_class_counts = []
        for read_class in classes_to_counts:
            read_classes.append(read_class)
            read_class_counts.append(classes_to_counts[read_class])
        reads_data = (read_classes, read_class_counts)
        self.output_miso_results(output_file,
                                 gene,
                                 reads_data,
                                 assignments,
                                 psi_vectors,
                                 kept_log_scores,
                                 num_iters,
                                 burn_in,
                                 lag,
                                 percent_acceptance,
                                 proposal_type)


    def output_miso_results(self, output_file, gene, reads_data, assignments,
                            psi_vectors, kept_log_scores, num_iters, burn_in,
                            lag, percent_acceptance, proposal_type):
        """
        Output results of MISO to a file.
        """
        output = open(output_file, 'w')

        # Get a string representation of the isoforms - use '_'
        # in the delimiter regardless
        iso_delim = '_'
        if type(gene.isoforms[0].desc) == list:
            str_isoforms = '[' + ",".join(["\'" + iso_delim.join(iso.desc) + "\'" \
                                           for iso in gene.isoforms]) + ']'
        else:
            str_isoforms = '[' + ",".join(["\'" + iso.desc + "\'" \
                                           for iso in gene.isoforms]) + ']'

        num_isoforms = len(gene.isoforms)

        # And of the exon lengths
        exon_lens = ",".join(["(\'%s\',%d)" %(p.label, p.len) \
                              for p in gene.parts])

        ## Compile header with information about isoforms and internal parameters used
        ## by the sampler, and also information about read counts and number of
        ## reads assigned to each isoform.
        read_classes, read_class_counts = reads_data
        read_counts_list = []

        for class_num, class_type in enumerate(read_classes):
            class_counts = read_class_counts[class_num]

            # Get the read class type in string format
            class_str = str(tuple([int(c) for c in class_type])).replace(" ", "")

            # Get the read class counts in string format
            class_counts_str = "%s" %(int(read_class_counts[class_num]))

            # Put class and counts together
            curr_str = "%s:%s" %(class_str,
                                 class_counts_str)
            read_counts_list.append(curr_str)

        # Get a summary of the raw read counts supporting each isoform
        read_counts_str = ",".join(read_counts_list)

        assigned_counts = reads_utils.count_isoform_assignments(assignments)

        # Get number of reads assigned to each isoform
        assigned_counts_str = ",".join(["%d:%d" %(c[0], c[1]) \
                                        for c in assigned_counts])

        # coordinates where mRNAs start
        mRNA_starts = []
        mRNA_ends = []
        for iso in gene.isoforms:
            mRNA_starts.append(iso.genomic_start)
            mRNA_ends.append(iso.genomic_end)
        mRNA_start_coords = ",".join([str(start) for start in mRNA_starts])
        mRNA_end_coords = ",".join([str(end) for end in mRNA_ends])
        chrom = gene.chrom
        if chrom == None:
            chrom = "NA"
        strand = gene.strand
        if strand == None:
            strand = "NA"
        header = "#isoforms=%s\texon_lens=%s\titers=%d\tburn_in=%d\tlag=%d\t" \
                 "percent_accept=%.2f\tproposal_type=%s\t" \
                 "counts=%s\tassigned_counts=%s\tchrom=%s\tstrand=%s\tmRNA_starts=%s\tmRNA_ends=%s\n" \
                 %(str_isoforms, exon_lens, num_iters, burn_in, lag,
                   percent_acceptance, proposal_type, read_counts_str,
                   assigned_counts_str,
                   # Fields related to gene/event
                   chrom,
                   strand,
                   mRNA_start_coords,
                   mRNA_end_coords)
        output.write(header)

        # Output samples and their associated log scores, as well as read counts
        results_fields = ["sampled_psi", "log_score"]
        results_header = "%s\n" %("\t".join(results_fields))
        output.write(results_header)
        for psi_sample, curr_log_score in zip(psi_vectors, kept_log_scores):
            psi_sample_str = ",".join(["%.4f" %(psi) for psi in psi_sample])
            output_line = "%s\t%.2f\n" %(psi_sample_str, curr_log_score)
            output.write(output_line)
        output.close()
        print "Completed outputting."
#        return [percent_acceptance, array(psi_vectors), array(kept_log_scores)]


def run_sampler_on_event(gene, ni, ne, nb, read_len, overhang_len, num_iters,
                         output_dir, confidence_level=.95):
    """
    Run sampler on a two-isoform gene event.
    """
    print "Running sampler on a two-isoform event..."
    print "  - Gene label: ", gene.label, gene
    print "  - NI, NE, NB: %d, %d, %d" %(ni, ne, nb)
    print "Using default sampler parameters."
    if gene.chrom != None:
        # Index output by chromosome
        print "Indexing by chromosome..."
        output_dir = os.path.join(output_dir, gene.chrom)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    output_filename = os.path.join(output_dir, gene.label)

    samples = []
    cred_interval = []

    num_isoforms = len(gene.isoforms)
    burn_in = 500
    lag = 10
    hyperparameters = ones(num_isoforms)
    proposal_diag = 0.05
    sigma = set_diag(zeros([num_isoforms-1, num_isoforms-1]),
                     proposal_diag)
    sampler_params = {'read_len': read_len,
                      'overhang_len': overhang_len,
                      'uniform_proposal': False,
                      'sigma_proposal': sigma}
    sampler = MISOSampler(sampler_params, log_dir=output_dir)
    reads = read_counts_to_read_list(ni, ne, nb)
    t1 = time.time()
    sampler_results = \
        sampler.run_sampler(num_iters, reads, gene, hyperparameters,
                            sampler_params, output_filename, burn_in=burn_in,
                            lag=lag)
    if not sampler_results:
        return (samples, cred_interval)
    samples = sampler_results[1]
    # Compute credible intervals
    cred_interval = ht.compute_credible_intervals(samples, confidence_level=confidence_level)
    t2 = time.time()
    print "  - Sampler run took %s seconds." %(str(t2-t1))
    # return samples and credible intervals
    return (samples, cred_interval)


def profile_miso():
    from Gene import make_gene
    gene = make_gene([150, 100, 150], [[1, 2, 3], [1, 3]])
    read_len = 36
    overhang_len = 4
    output_dir = "profiler-test"
    for x in range(10):
        print "x = %d" %(x)
        a, b = run_sampler_on_event(gene, 500, 50, 40, read_len, overhang_len,
                                    10000, output_dir)



def main():
    return
    # import cProfile as profile
    # import pstats
    # output_file = "profile"
    # profile.run('profile_miso()', output_file)
    # p = pstats.Stats(output_file)
    # print "name: "
    # print p.sort_stats('name')
    # print "all stats: "
    # p.print_stats()
    # print "cumulative (top 10): "
    # p.sort_stats('cumulative').print_stats(20)

if __name__ == '__main__':
    main()
