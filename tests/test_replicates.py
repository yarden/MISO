
import setup_miso
import misopy
import misopy.gff_utils
import pysplicing
import random as plainrandom

from scipy import *
from numpy import *

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def do_replicates(reads):
    result = pysplicing.doMISO(
        GFF = gene, gene = 0L, reads = reads, read_len = 33L,
        num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC,
        replicate_mean_prior_mean = 0.0, replicate_mean_prior_var = 100.0,
        replicate_var_prior_numobs = 0.01, replicate_var_prior_var = 0.01)
    psi = transpose(array(result[0]))
    meanpsi = mean(psi, 0)
    return meanpsi

def do_one(reads):
    results = pysplicing.doMISO(
        gene, 0L, (reads,), read_len = 33L, num_iters = 5000L,
        burn_in = 1000L, lag = 10L)
    psi = transpose(array(results[0]))
    meanpsi = mean(psi, 0)
    return meanpsi

def do_onebyone(reads):
    psi = matrix([ do_one(r) for r in reads ])
    return mean(psi, 0).ravel().tolist()[0]

def test_replicates():
    plainrandom.seed(42)
    reads1 = pysplicing.simulateReads(gene, 0L, (0.1, 0.9), 100L, 33L)
    reads2 = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 100L, 33L)
    reads3 = pysplicing.simulateReads(gene, 0L, (0.3, 0.7), 100L, 33L)
    reads = (tuple(reads1[1:3]), tuple(reads2[1:3]), tuple(reads3[1:3]))
    psi_rep = do_replicates(reads)
    psi_one = do_onebyone(reads)

    print psi_rep
    print psi_one
