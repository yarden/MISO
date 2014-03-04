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
import misopy.sam_utils as sam_utils
import misopy.reads_utils as reads_utils

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


    def test_a_gene_class(self):
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
        ##
        ## Test gene loading directly from GFF
        ##
        print "Loading gene from GFF"
        gene_obj.from_gff(gff_fname)
        # Now load it with a name parameter
        gene_obj.from_gff(gff_fname, gene_id="ENSMUSG00000019943")
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
          "Cannot retrieve dicitonary from JSON"
        # Create new gene object with same properties
        new_gene_obj = gene_class.Gene(from_json_fname="./__gene.json")
        print "New gene made from JSON"
        assert (new_gene_obj == gene_obj), "New gene copy not equal to old."


    def test_b_alignments(self):
        """
        Test alignments of reads to genes.
        """
        # Example of an SE event from mm9
        example = test_cases.get_example("mm9_se_example")
        # First load the GFF as a Gene object
        gff_fname = example["gff"]
        gene_obj = gene_class.Gene(from_gff_fname=gff_fname)
        print "Gene obj: ", gene_obj
        bam_fname = example["bam"]
        # Load BAM file
        bamfile = sam_utils.load_bam_reads(bam_fname)
        # Fetch reads in gene
        bam_reads = sam_utils.fetch_bam_reads_in_gene(bamfile, gene_obj.chrom,
                                                      gene_obj.start,
                                                      gene_obj.end)
        read_len = 48
        overhang_len = 1
        alignments = gene_obj.align_single_end_reads(bam_reads, read_len,
                                                     overhang_len=overhang_len)
        # Expected counts for 'sample1' are: counts=(0,1):1,(1,0):21,(1,1):23
        expected_counts = {(0, 1): 1,
                           (1, 0): 21,
                           (1, 1): 23}
        read_class_counts = reads_utils.count_aligned_reads(alignments)
        read_class_counts = dict(read_class_counts)
        for counts_class in expected_counts:
            assert (read_class_counts[counts_class] == \
                    expected_counts[counts_class]), \
                    "Failed to get right number for %s. Got %d, expected " \
                    "%d." %(counts_class,
                            read_class_counts[counts_class],
                            expected_counts[counts_class])
        
        


    def test_pairing(self):
        """
        Test pairing of reads.
        """
        pass


def main():
    unittest.main()


if __name__ == "__main__":
    main()
