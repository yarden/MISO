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

import misopy.Gene as gene_utils


class TestGene(unittest.TestCase):
    """
    Test gene.
    """
    def setUp(self):
        pass


    def test_gene_class(self):
        label = "mygene"
        chrom = "chr22"
        strand = "+"
        isoform_desc = [["A", "B", "C"], ["A", "C"]]
        # Set gene programmatically
        exons = [gene_utils.Exon(10,20), gene_utils.Exon(30,40), gene_utils.Exon(100,200)]
        gene_obj = gene_class.Gene(label, chrom, strand)
        print "Gene obj: ", gene_obj
        # print gene_obj.label
        # print gene_obj.isoform_desc
        # for iso in gene_obj.isoform_desc:
        #     print "Isoform: ", iso
        # for exon in gene_obj.parts:
        #     print "Exon: ", exon
        # # Overwrite exons
        # new_exons = [gene_utils.Exon(5,20)]
        # gene_obj.parts = new_exons
        # print "Resetting exons: "
        # print gene_obj.parts



def main():
    unittest.main()


if __name__ == "__main__":
    main()

