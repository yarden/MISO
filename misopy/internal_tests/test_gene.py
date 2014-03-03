##
## Test gene utilities
## 
## Yarden Katz <yarden@mit.edu>
##
import os
import sys
import time
import unittest

import numpy as np
import numpy.linalg as linalg
import math

import misopy
import misopy.internal_tests
import misopy.internal_tests.py_scores as py_scores

import misopy.pyx
import misopy.pyx.miso_scores_single as scores_single
import misopy.pyx.miso_scores_paired as scores_paired
import misopy.pyx.stat_helpers as stat_helpers
import misopy.pyx.math_utils as math_utils
import misopy.pyx.gene_class as gene_class


class TestScores(unittest.TestCase):
    """
    Test MISO scoring functions.
    """
    def setUp(self):
        pass


    def test_gene_class(self):
        label = "mygene"
        chrom = "chr22"
        strand = "+"
        isoform_desc = [["A", "B", "C"], ["A", "C"]]
        exons = [["A", 100], ["B", 50], ["C", 120]]
        gene_obj = gene_class.MISOGene(label, chrom, strand, isoform_desc, exons)
        print "Gene obj: ", gene_obj
        print gene_obj.label
        print gene_obj.isoform_desc
        for iso in gene_obj.isoform_desc:
            print "Isoform: ", iso
        for exon in gene_obj.exons:
            print "Exon: ", exon
            print exon[0], exon[1]



def main():
    unittest.main()


if __name__ == "__main__":
    main()

