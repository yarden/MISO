
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
    random.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)

    random.seed(42)
    results1 = pysplicing.doMISO(
        gene, 0L, reads[1], reads[2], 33L, 500L, 100L, 10L,
        dirichlet_prior_params = (1.0, 1.0))
    psi1 = transpose(array(results1[0]))

    random.seed(42)
    results2 = pysplicing.doMISO(
        GFF = gene, gene = 0L, read_pos = reads[1], read_cigar = reads[2],
        read_len = 33L, num_iters = 500L, burn_in = 100L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC, logistic_prior_mean = 0.0,
        logistic_prior_var = 3.0)
    psi2 = transpose(array(results2[0]))
