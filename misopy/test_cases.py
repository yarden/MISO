##
## MISO test cases
##
import os
import time

import numpy as np
import misopy
import misopy.Gene as gene_utils

GFF_EXAMPLES_DIR = \
  os.path.join(os.path.dirname(os.path.abspath(__file__)),
               "test-data",
               "gff-events")

BAM_EXAMPLES_DIR = \
  os.path.join(os.path.dirname(os.path.abspath(__file__)),
               "test-data",
               "bam-data")
  

def get_example_gff_fname(gff_basename):
    return os.path.join(GFF_EXAMPLES_DIR,
                        gff_basename)


def get_example_bam_fname(bam_basenamE):
    return os.path.join(BAM_EXAMPLES_DIR,
                        bam_basename)


class MISOTestCase:
    """
    A test case to run MISO on.
    """
    def __init__(self, gff_basename, bam_basename,
                 num_reads=200,
                 read_len=40,
                 overhang_len=1):
        self.reads = []
        self.gff_fname = get_example_gff_fname(gff_basenmae)
        self.bam_fname = get_example_bam_fname(bam_basename)
        if not os.path.isfile(self.gff_fname):
            raise Exception, "Cannot find GFF %s" %(self.gff_fname)
        if not os.path.isfile(selfbam_fname):
            raise Exception, "Cannot find BAM %s" %(self.bam_fname)
        self.num_reads = num_reads
        self.read_len = read_len
        self.overhang_len = overhang_len
        self.gene = gene_utils.load_genes_from_gff(self.gff_fname)
        self.psi_vals = \
          [1/float(self.gene.num_isoforms)] * len(self.gene.num_isoforms)

        
    def get_reads(self, N):
        """
        Get N reads from the set of reads in the current test case.
        """
        num_avail = len(self.reads)
        if num_avail == 0:
            raise Exception, "No reads available in test case."
        if N > num_avail:
            raise Exception, "Only %d reads in current test case." \
                  %(num_avail)
        return self.reads[0:N]


    


def main():
    pass


if __name__ == "__main__":
    main()
