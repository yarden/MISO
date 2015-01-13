
import misopy
import misopy.gff_utils
# TODO: this loads the _installed_ pysplicing
import pysplicing
import random

from scipy import *
from numpy import *

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def test_MISO():
    reads = pysplicing.simulateReads(gene, 0L, (0.5, 0.5), 10L, 33L)

    results1 = pysplicing.doMISO(
        gene, 0L, reads[1], reads[2], 33L, 5000L, 1000L, 10L)
    psi1 = transpose(array(results1[0]))

    print
    print "Dirichlet:", mean(psi1, 0)

    results2 = pysplicing.doMISO(
        GFF = gene, gene = 0L, read_pos = reads[1], read_cigar = reads[2],
        read_len = 33L, num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = 0.0,
        logistic_prior_var = .01)
    psi2 = transpose(array(results2[0]))

    print "Logistic: ", mean(psi2, 0)

    results3 = pysplicing.doMISO(
        GFF = gene, gene = 0L, read_pos = reads[1], read_cigar = reads[2],
        read_len = 33L, num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = 5.0,
        logistic_prior_var = 1.0)
    psi3 = transpose(array(results3[0]))

    print "Logistic: ", mean(psi3, 0)

if __name__ == "__main__":
    test_MISO()
