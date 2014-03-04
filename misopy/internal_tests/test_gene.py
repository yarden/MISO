##
## Test gene utilities
## 
## Yarden Katz <yarden@mit.edu>
##
import os
import sys
import time
import unittest
import copy

import numpy as np
import numpy.linalg as linalg
import math

import misopy
import misopy.internal_tests
import misopy.internal_tests.py_scores as py_scores
import misopy.internal_tests.test_cases as test_cases

import misopy.gff_utils as gff_utils

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
        # Load the GFF here
        self.test_data_dir = test_cases.get_test_data_dir()


    def test_gene_class(self):
        label = "mygene"
        chrom = "chr22"
        strand = "+"
        isoform_desc = [["A", "B", "C"], ["A", "C"]]
        # Set gene programmatically
        exons = [gene_utils.Exon(10,20), gene_utils.Exon(30,40), gene_utils.Exon(100,200)]
        gene_obj = gene_class.Gene(label, chrom, strand, 10, 20)
        print "Gene obj: ", gene_obj
        # Make gene from GFF
        gff_fname = test_cases.get_gene_gff("Atp2b1")
        gff_db = gff_utils.GFFDatabase(gff_fname)
        gene_id = gff_db.genes[0].get_id()
        results = gff_db.get_genes_records([gene_id])
        gene_recs = results[0]
        gene_hierarchy = results[1][gene_id]
        print "Loading gene from GFF"
        gene_obj.from_gff_recs(gene_hierarchy, gene_recs)
        # Take first part of first transcript
        part = gene_obj.transcripts[0].parts[0]
        # Check that first transcript has that part
        assert (gene_obj.transcripts[0].has_part(part)), \
          "Cannot find part in transcript."
        # Make a copy and check for presence of part
        part_copy = copy.copy(part)
        assert (gene_obj.transcripts[0].has_part(part_copy)), \
          "Cannot find part copy in transcript."
        for transcript in gene_obj.transcripts:
            print "Transcript:" , transcript.label
        # Test retrieval of constitutive parts
        const_parts = gene_obj.get_const_parts()
        assert (len(const_parts) == 10), "Failed to get all constitutive exons."
        ##
        ## Test gene serialization and retrieval with JSON
        ##
        # Serialize as JSON
        gene_obj.save_json("./__gene.json")
        # Read gene back from JSON
        new_gene_dict = gene_obj.from_json("./__gene.json")
        assert (new_gene_dict is not None), \
          "Cannot retriev dicitonary from JSON"
        # Create new gene object with same properties
        new_gene_obj = gene_class.Gene(from_json_fname="./__gene.json")
        print "New gene made from JSON"
        assert (new_gene_obj == gene_obj), "New gene copy not equal to old."


def main():
    unittest.main()


if __name__ == "__main__":
    main()
