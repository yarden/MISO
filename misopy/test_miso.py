#!/usr/bin/env python
import os
import sys
import unittest

import pysam
import sam_utils

import misopy
import misopy.gff_utils as gff_utils
import misopy.reads_utils as reads_utils

class TestMISO(unittest.TestCase):
    """
    Test MISO functionality.
    """
    def setUp(self):
        # Find out the current directory
        self.miso_path = \
            os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
        self.tests_data_dir = \
            os.path.join(self.miso_path, "test-data")
        self.events_analysis_cmd = \
            "python %s " %(os.path.join(self.miso_path,
                                        "run_events_analysis.py"))
        self.tests_output_dir = \
            os.path.join(self.miso_path, "test-output")
        self.test_sam_filename = \
            os.path.join(self.tests_data_dir,
                         "sam-data",
                         "c2c12.Atp2b1.sam")
        self.gff_events_dir = \
            os.path.join(self.miso_path, "gff-events")
        self.sam_to_bam_script = \
            os.path.join(self.miso_path, "sam_to_bam.py")
        self.index_gff_script = \
            os.path.join(self.miso_path, "index_gff.py")
        # Test event BAM file and its index
        self.test_event_bam_fname = \
            os.path.join(self.tests_data_dir, "sam-data", "test_event_reads.bam")
        self.test_event_bam_index_fname = \
            "%s.bai" %(os.path.join(self.test_event_bam_fname))
        


    def ___test_a_sam_to_bam(self):
        """
        Test conversion of SAM to BAM.

        The 'a' ensures this runs first.
        """

        print "Testing conversion of SAM to BAM..."
        output_dir = \
            os.path.join(self.tests_output_dir, "sam-output")
        sam_to_bam_cmd = \
            "python %s --convert %s %s" %(self.sam_to_bam_script,
                                          self.test_sam_filename,
                                          output_dir)
        print "Executing: %s" %(sam_to_bam_cmd)
        os.system(sam_to_bam_cmd)

        # Make sure conversion worked; sorted, indexed BAM file is outputted
        assert(os.path.exists(os.path.join(output_dir,
                                           "c2c12.Atp2b1.sorted.bam")))


    def test_a2_strandedness(self):
        """
        Test that strandedness is read correctly.
        """
        # Read 1 is forward, on plus strand
        # Has flag 129, i.e. '0b10000001'
        f_read = pysam.AlignedRead()
        f_read.qname = "f_read"
        f_read.flag = 129
        f_read.rname = 9
        f_read.pos = 4991443

        # Read 2 is reverse, on minus strand
        # Has flag 81, i.e. '0b1010001'
        r_read = pysam.AlignedRead()        
        r_read.qname = "r_read"
        r_read.flag = 81
        r_read.rname = 9
        r_read.pos = 4991578

        # Test that we can read the BAM strand flag correctly
        assert(sam_utils.flag_to_strand(f_read.flag) == "+"), \
            "Error in determining plus strand of read."
        assert(sam_utils.flag_to_strand(r_read.flag) == "-"), \
            "Error in determining minus strand of read."
        ##
        ## Test stranded-ness rules
        ##
        #   fr-unstranded,
        #   fr-firststrand,
        #   fr-secondstrand
        plus_target_strand = "+"
        minus_target_strand = "-"
        # fr-unstranded: both strand reads should match
        # either target strand
        print "Testing fr-unstranded..."
        for curr_read in [f_read, r_read]:
            for target in [plus_target_strand, minus_target_strand]:
                print "Checking read ", curr_read.qname, " against ", target
                assert(sam_utils.read_matches_strand(curr_read,
                                                     target,
                                                     "fr-unstranded") == True), \
                    "Error checking strand of fr-unstranded."
        # fr-firststrand: forward read must be *opposite* of target strand,
        # i.e. +read matches -target, and -read matches +target
        # test +read
        print "Testing fr-firststrand..."
        assert(sam_utils.read_matches_strand(f_read,
                                             minus_target_strand,
                                             "fr-firststrand") == True), \
            "+read must match -target under fr-firstrand."
        assert(sam_utils.read_matches_strand(f_read,
                                             plus_target_strand,
                                             "fr-firststrand") == False), \
            "+read must *not* match +target under fr-firststrand."
        # test -read
        assert(sam_utils.read_matches_strand(r_read,
                                             plus_target_strand,
                                             "fr-firststrand") == True), \
            "-read must match +target under fr-firststrand."
        assert(sam_utils.read_matches_strand(r_read,
                                             minus_target_strand,
                                             "fr-firststrand") == False), \
            "+read must match -target under fr-firststrand."
        # Test fr-firststrand read pair
        pe = (300, 10)
        assert(sam_utils.read_matches_strand((f_read, r_read),
                                             minus_target_strand,
                                             "fr-firststrand",
                                             paired_end=pe) == True), \
            "(+, -) must match -target under fr-firststrand."
        assert(sam_utils.read_matches_strand((f_read, r_read),
                                             plus_target_strand,
                                             "fr-firststrand",
                                             paired_end=pe) == False), \
            "(+, -) must *not* match +target under fr-firststrand."
        # fr-secondstrand: forward read must be *matching* target strand,
        # i.e. +read matches +target, and -read matches +target
        print "Testing fr-secondstrand..."
        # test +read
        assert(sam_utils.read_matches_strand(f_read,
                                             plus_target_strand,
                                             "fr-secondstrand") == True), \
            "+read must match +target under fr-secondstrand."
        assert(sam_utils.read_matches_strand(f_read,
                                             minus_target_strand,
                                             "fr-secondstrand") == False), \
            "+read must *not* match -target under fr-secondstrand."
        # test -read
        assert(sam_utils.read_matches_strand(r_read,
                                             minus_target_strand,
                                             "fr-secondstrand") == True), \
            "-read must match -target under fr-secondstrand."
        assert(sam_utils.read_matches_strand(r_read,
                                             plus_target_strand,
                                             "fr-secondstrand") == False), \
            "-read must *not* match +target under fr-secondstrand."
        # Test fr-secondstrand read pair
        assert(sam_utils.read_matches_strand((f_read, r_read),
                                             plus_target_strand,
                                             "fr-secondstrand",
                                             paired_end=pe) == True), \
            "(+, -) must match + under fr-secondstrand."
        assert(sam_utils.read_matches_strand((f_read, r_read),
                                             minus_target_strand,
                                             "fr-secondstrand",
                                             paired_end=pe) == False), \
            "(+, -) must *not* match - under fr-secondstrand."                                 

        
    def ___test_z_gene_psi(self):
        """
        Test gene-level Psi inferences using SAM/BAM reads.

        The 'z' ensures this runs last.
        """
        print "Testing gene-level Psi..."
        sam_dir = os.path.join(self.tests_output_dir, "sam-output")
        bam_filename = os.path.join(sam_dir, "c2c12.Atp2b1.sorted.bam")

        read_len = 35
        insert_mean = 250
        insert_sd = 30

        # First index the GFF of interest
        gff_filename = os.path.join(self.gff_events_dir,
                                    "mm9",
                                    "genes",
                                    "Atp2b1.mm9.gff")
        gff_index_dir = os.path.join(self.gff_events_dir,
                                     "mm9",
                                     "genes",
                                     "Atp2b1",
                                     "indexed")
        print "Testing GFF indexing of: %s" %(gff_filename)
        index_cmd = "python %s --index %s %s" %(self.index_gff_script,
                                                gff_filename,
                                                gff_index_dir)

        print "Executing: %s" %(index_cmd)
        os.system(index_cmd)

        output_dir = os.path.join(self.tests_output_dir,
                                  "gene-psi-output")
        miso_cmd = "%s --compute-genes-psi %s %s --output-dir %s --read-len %d " \
                   %(self.events_analysis_cmd,
                     gff_index_dir,
                     bam_filename,
                     output_dir,
                     read_len)
        print "Executing: %s" %(miso_cmd)
        os.system(miso_cmd)


    def test_isoform_assignments(self):
        """
        Test isoform assignments.
        """
        print "Testing isoform assignments..."
        gff_filename = os.path.join(self.gff_events_dir,
                                    "mm9",
                                    "events",
                                    "test_event.gff")
        gff_index_dir = os.path.join(self.gff_events_dir,
                                     "mm9",
                                     "events",
                                     "test_event_indexed")
        print "Making indexed GFF from: %s" %(gff_filename)
        index_cmd = "python %s --index %s %s" %(self.index_gff_script,
                                                gff_filename,
                                                gff_index_dir)
        # Index the event
        os.system(index_cmd)
        # BAM file to run on: check that it exists
        if not os.path.isfile(self.test_event_bam_fname):
            raise Exception, "Cannot find test event BAM file %s" \
                  %(self.test_event_bam_fname)
        # Check that the index is there too
        if not os.path.isfile(self.test_event_bam_index_fname):
            raise Exception, "Cannot find test event BAM index file %s" \
                  %(self.test_event_bam_index_fname)
        # Load the gene from the GFF
        test_event_chrom = \
            "chr11"
        test_event_id = \
            "chr11:77667254:77667393:+@chr11:77668180:77668224:+@chr11:77671225:77671347:+"
        gff_index_filename = \
            os.path.join(gff_index_dir, "chr11", "%s.pickle" %(test_event_id))
        if not os.path.isfile(gff_index_filename):
            raise Exception, "Cannot find pickle file for test_event.gff."
        gff_genes = gff_utils.load_indexed_gff_file(gff_index_filename)
        # Retrieve gene
        gene_id = gff_genes.keys()[0]
        gene_info = gff_genes[gene_id]
        gene_obj = gene_info['gene_object']
        gene_hierarchy = gene_info['hierarchy']
        # Find the most inclusive transcription start and end site
        # for the gene
        print "Retrieving gene start and end boundaries..."
        tx_start, tx_end = \
            gff_utils.get_inclusive_txn_bounds(gene_info['hierarchy'][gene_id])
        # Fetch reads aligning to the gene boundaries
        print "Fetching BAM reads in gene..."
        bamfile = pysam.Samfile(self.test_event_bam_fname, "rb")
        gene_reads = \
            sam_utils.fetch_bam_reads_in_gene(bamfile,
                                              gene_obj.chrom,
                                              tx_start,
                                              tx_end,
                                              gene_obj)
        print "Fetched reads successfully."
        # Align reads to isoform
        read_len = 40
        overhang_len = 1
        print "Aligning reads to isoforms..."
        if len(gene_obj.isoforms) != 2:
            raise Exception, "Unexpected number of isoforms (needed 2.)"
        read_assignments = \
            sam_utils.sam_reads_to_isoforms(gene_reads, gene_obj, read_len,
                                            overhang_len)
        print "Counting number of reads in each read class..."
        print "reads to isoforms: "
        print read_assignments
        # Count isoforms assignments
        print "Counting assignments of reads to isoforms..."
        ##
        ## Correct counts should be:
        ##
        # counts=(0,0):222,(0,1):21,(1,0):63,(1,1):272        
        counts = reads_utils.collapse_isoform_assignments(read_assignments)
        print "Counts: "
        print counts

        
if __name__ == '__main__':
    unittest.main()
