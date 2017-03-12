
import setup_miso
import misopy
import misopy.gff_utils
import pysplicing
import random as plainrandom

from scipy import *
from numpy import *

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def test_MISO():
    plainrandom.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.5, 0.5), 10L, 33L)

    results1 = pysplicing.doMISO(
        gene, 0L, (tuple(reads[1:3]),), 33L, 5000L, 1000L, 10L)
    psi1 = transpose(array(results1[0]))
    meanpsi1 = mean(psi1, 0)

    print "Dirichlet:", meanpsi1

    assert( abs(meanpsi1[0] - 0.43982318) < 1e-8 and
            abs(meanpsi1[1] - 0.56017682) < 1e-8 )

    results2 = pysplicing.doMISO(
        GFF = gene, gene = 0L, reads = (tuple(reads[1:3]),),
        read_len = 33L, num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = 0.0,
        logistic_prior_var = .01)
    psi2 = transpose(array(results2[0]))
    meanpsi2 = mean(psi2, 0)

    print "Logistic: ", meanpsi2

    assert( abs(meanpsi2[0] - 0.49930664) < 1e-8 and
            abs(meanpsi2[1] - 0.50069336) < 1e-8 )

    results3 = pysplicing.doMISO(
        GFF = gene, gene = 0L, reads = (tuple(reads[1:3]),),
        read_len = 33L, num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = 5.0,
        logistic_prior_var = 1.0)
    psi3 = transpose(array(results3[0]))
    meanpsi3 = mean(psi3, 0)

    print "Logistic: ", meanpsi3

    assert( abs(meanpsi3[0] - 0.92413627) < 1e-8 and
            abs(meanpsi3[1] - 0.07586373) < 1e-8 )

    results4 = pysplicing.doMISO(
        GFF = gene, gene = 0L, reads = (tuple(reads[1:3]),),
        read_len = 33L, num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = -5.0,
        logistic_prior_var = 1.0)
    psi4 = transpose(array(results4[0]))
    meanpsi4 = mean(psi4, 0)

    print "Logistic: ", meanpsi4

    assert( abs(meanpsi4[0] - 0.08581447) < 1e-8 and
            abs(meanpsi4[1] - 0.91418553) < 1e-8 )

def test_paired_MISO():
    plainrandom.seed(42)
    reads = pysplicing.simulatePairedReads(gene,       # gff
                                           0L,         # gene
                                           (0.5, 0.5), # expression
                                           10L,        # noreads
                                           33L,        # readlength
                                           60,         # normalMean
                                           5,          # normalVar
                                           4)          # numDevs

    results1 = pysplicing.doMISOPaired(
        gene, 0L, (tuple(reads[1:3]),), 33L, 60, 5, 4, 5000L, 1000L, 10L)
    psi1 = transpose(array(results1[0]))
    meanpsi1 = mean(psi1, 0)

    print "Dirichlet:", meanpsi1

    assert( abs(meanpsi1[0] - 0.36982598) < 1e-8 and
            abs(meanpsi1[1] - 0.63017402) < 1e-8 )

    results2 = pysplicing.doMISOPaired(
        GFF = gene, gene = 0L, reads = (tuple(reads[1:3]),),
        read_len = 33L, mean_frag_len = 60, frag_variance = 5, num_sds = 4,
        num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = 0.0,
        logistic_prior_var = 0.01)
    psi2 = transpose(array(results2[0]))
    meanpsi2 = mean(psi2, 0)

    print "Logistic: ", meanpsi2

    assert( abs(meanpsi2[0] - 0.39687208) < 1e-8 and
            abs(meanpsi2[1] - 0.60312792) < 1e-8 )

    results3 = pysplicing.doMISOPaired(
        GFF = gene, gene = 0L, reads = (tuple(reads[1:3]),),
        read_len = 33L, mean_frag_len = 60, frag_variance = 5, num_sds = 4,
        num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = 5.0,
        logistic_prior_var = 1.0)
    psi3 = transpose(array(results3[0]))
    meanpsi3 = mean(psi3, 0)

    print "Logistic: ", meanpsi3

    assert( abs(meanpsi3[0] - 0.39861373) < 1e-8 and
            abs(meanpsi3[1] - 0.60138627) < 1e-8 )

if __name__ == "__main__":
    test_MISO()
    test_paired_MISO()
