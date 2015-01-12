
import misopy
import misopy.gff_utils
# TODO: this loads the _installed_ pysplicing
import pysplicing
import random

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def test_rng_simulate_reads():
    random.seed(42)
    reads1 = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    random.seed(42)
    reads2 = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    assert(reads1 == reads2)

def test_rng_MISO():
    random.seed(42)
    reads1 = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    results1 = pysplicing.MISO(gene, 0L, reads1[1], reads1[2], 33L, 50L,
                               10L, 2L, (1.0, 1.0))
    random.seed(42)
    reads2 = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    results2 = pysplicing.MISO(gene, 0L, reads2[1], reads2[2], 33L, 50L,
                               10L, 2L, (1.0, 1.0))
    assert(results1 == results2)
