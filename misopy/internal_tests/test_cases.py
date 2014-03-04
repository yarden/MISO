##
## Test cases
##
import os
import sys
import time

import numpy as np

import misopy

num_inc = 3245
num_exc = 22
num_com = 39874
READS = [[1,0]] * num_inc + \
        [[0,1]] * num_exc + \
        [[1,1]] * num_com
READS = np.array(READS, dtype=np.dtype("i"))
read_len = 40
overhang_len = 4
num_parts_per_iso = np.array([3, 2], dtype=np.dtype("i"))
iso_lens = np.array([1253, 1172], dtype=np.dtype("i"))
# Assignment of reads to isoforms: assign half of
# the common reads to isoform 0, half to isoform 1
iso_nums = [0]*3245 + [1]*22 + [0]*19937 + [1]*19937
iso_nums = np.array(iso_nums, dtype=np.dtype("i"))
num_reads = len(READS)


def get_test_data_dir():
    test_data_dir = os.path.abspath(os.path.join(__file__, "..",
                                                 "..", "test-data"))
    return test_data_dir


def get_bam_dir()
    data_dir = get_test_data_dir()
    return os.path.join(data_dir, "bam-data")


def get_gene_gff(example_name):
    print "Getting gene GFF for %s" %(example_name)
    data_dir = get_test_data_dir()
    fname = None
    if example_name == "Atp2b1":
        fname = os.path.join(data_dir, "gff-events", "mm9", 
                             "genes", "%s.mm9.gff" %(example_name))
    else:
        raise Exception, "Cannot find example %s" %(example_name)
    return fname


def get_bam_example(example_name):
    """
    Return an example of a gene/event GFF and an associated
    BAM file.
    """
    bam_dir = get_bam_dir()
    
    




