
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
        replicate_mean_prior_mean = 0.0, replicate_mean_prior_var = 100000,
        replicate_var_prior_numobs = 0.0001, replicate_var_prior_var = 0.0001,
        algorithm = pysplicing.MISO_ALGO_REASSIGN)
    psi = transpose(array(result[0]))
    meanpsi = mean(psi, 0)
    varpsi = var(psi, 0)
    assert(abs(varpsi[0] - 0.00696644) < 1e-8 and
           abs(varpsi[1] - 0.00696644) < 1e-8)
    return meanpsi

def do_one(reads):
    results = pysplicing.doMISO(
        gene, 0L, (reads,), read_len = 33L, num_iters = 5000L,
        burn_in = 1000L, lag = 10L, algorithm = pysplicing.MISO_ALGO_REASSIGN)
    psi = transpose(array(results[0]))
    return psi

def do_onebyone(reads):
    psi = [ do_one(r) for r in reads ]
    meanpsi = matrix([ mean(p, 0) for p in psi ])
    assert(var(meanpsi, 0).ravel().tolist()[0] ==
           [0.017807353098872165, 0.01780735309887213])
    return mean(meanpsi, 0).ravel().tolist()[0]

def test_replicates():
    plainrandom.seed(42)
    reads1 = pysplicing.simulateReads(gene, 0L, (0.1, 0.9), 100L, 33L)
    reads2 = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 100L, 33L)
    reads3 = pysplicing.simulateReads(gene, 0L, (0.3, 0.7), 100L, 33L)
    reads4 = pysplicing.simulateReads(gene, 0L, (0.4, 0.6), 100L, 33L)
    reads5 = pysplicing.simulateReads(gene, 0L, (0.5, 0.5), 100L, 33L)

    reads = (tuple(reads1[1:3]), tuple(reads2[1:3]), tuple(reads3[1:3]),
             tuple(reads4[1:3]), tuple(reads5[1:3]))
    psi_rep = do_replicates(reads)
    psi_one = do_onebyone(reads)

    assert(abs(psi_rep[0] - 0.31053177) < 1e-8 and
           abs(psi_rep[1] - 0.68946823) < 1e-8)

    assert(abs(psi_one[0] - 0.3232436493675063) < 1e-8 and
           abs(psi_one[1] - 0.6767563506324935) < 1e-8)
