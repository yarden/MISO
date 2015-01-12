
import misopy
import misopy.gff_utils
# TODO: this loads the _installed_ pysplicing
import pysplicing
import random

from scipy import *
from numpy import *

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def test_noIso():
    assert(pysplicing.noIso(gene) == (2,))

def test_isoLength():
    assert(pysplicing.isoLength(gene) == [[300, 200]])

def test_simulate_reads():
    random.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    assert(mean(reads[0]) == 0.7185)

def test_MISO():
    random.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    random.seed(42)
    results = pysplicing.MISO(gene, 0L, reads[1], reads[2], 33L, 500L,
                              100L, 10L)

    psi = transpose(array(results[0]))
    assert(all(mean(psi, 0) - [ 0.19308989, 0.80691011] < 1e-8))
