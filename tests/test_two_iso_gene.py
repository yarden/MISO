
import setup_miso
import misopy
import misopy.gff_utils
import pysplicing
import random as plainrandom

from scipy import *
from numpy import *

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def test_noIso():
    assert(pysplicing.noIso(gene) == (2,))

def test_isoLength():
    assert(pysplicing.isoLength(gene) == [[300, 200]])

def test_simulate_reads():
    plainrandom.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    meanreads = mean(reads[0])
    print(meanreads)
    assert(abs(meanreads - 0.7255) < 1e-8)

def test_MISO():
    def absl(x):
        return [ abs(e) for e in x ]

    plainrandom.seed(42)
    reads = pysplicing.simulateReads(gene, 0L, (0.2, 0.8), 2000L, 33L)
    results = pysplicing.MISO(gene, 0L, (tuple(reads[1:3]),), 33L, 500L,
                              100L, 10L)

    psi = transpose(array(results[0]))
    meanpsi = mean(psi, 0)
    diff = absl(meanpsi - [0.18950034, 0.81049966])
    print(meanpsi)
    print(diff)

    assert(diff[0] < 1e-8 and diff[1] < 1e-8)
