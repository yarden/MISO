
import misopy
import misopy.gff_utils
# TODO: this loads the _installed_ pysplicing
import pysplicing
import random

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def test_noIso():
    assert(pysplicing.noIso(gene) == (2,))

def test_isoLength():
    assert(pysplicing.isoLength(gene) == [[300, 200]])

def test_simulate_reads():
    random.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)

def test_MISO():
    random.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    random.seed(42)
    results = pysplicing.MISO(gene, 0L, reads[1], reads[2], 33L, 50L,
                              10L, 2L, (1.0, 1.0))
