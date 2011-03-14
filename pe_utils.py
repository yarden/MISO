##
## Utilities for working with paired-end reads and
## fragment distributions
##

import os
import glob
import time

from scipy import *
from numpy import *

from parse_csv import *
import Gene as gene_utils
import sam_utils


def load_insert_dist(insert_dist_filename):
    """
    Read insert length distribution.
    """
    insert_dist_file = open(insert_dist_filename, "r")
    insert_dist = array([int(line.strip()) \
                         for line in insert_dist_file])
    return insert_dist

    # insert_filename = glob.glob(os.path.join(insert_dir, "*.insert_dist"))
    # if len(insert_filename) == 0:
    #     raise Exception, "No insert length file found in %s." %(insert_dir)


def compute_insert_len(bam_filename, gff_filename, output_dir,
                       min_exon_size):
    """
    Compute insert length distribution and output it to the given
    directory.
    """
    print "Computing insert length distribution of %s" %(bam_filename)
    print "  - Using gene models from: %s" %(gff_filename)
    print "  - Outputting to: %s" %(output_dir)
    print "  - Minimum exon size used: %d" %(min_exon_size)

    if not os.path.isdir(output_dir):
        print "Making directory: %s" %(output_dir)
        os.makedirs(output_dir)

    output_filename = os.path.join(output_dir,
                                   "%s.insert_len" %(os.path.basename(bam_filename)))

    # Load BAM file with reads
    bamfile = sam_utils.load_bam_reads(bam_filename)
    
    # Load the genes from the GFF
    print "Loading genes from GFF..."
    t1 = time.time()
    gff_genes = gene_utils.load_genes_from_gff(gff_filename)
    t2 = time.time()
    print "  - Loading genes from GFF took %.2f seconds" %(t2 - t1)

    insert_lengths = []

    t1 = time.time()

    relevant_region = 0
    
    for gene_id, gene_info in gff_genes.iteritems():
        gene_obj = gene_info["gene_object"]

        # Get all the constitutive parts
        const_parts = gene_obj.get_const_parts()

        chrom = gene_obj.chrom

        # Consider only the large constitutive parts
        for part in const_parts:
            if part.len >= min_exon_size:
                # Get all the reads that land in the coordinates of the exon
                try:
                    exon_reads = bamfile.fetch(chrom, part.start, part.end)
                except ValueError:
                    print "Could not fetch from region: ", chrom, part.start, part.end
                    continue

                # Pair all the paired-end reads that land there
                paired_reads = sam_utils.pair_sam_reads(exon_reads)
                num_paired_reads = len(paired_reads)

                if num_paired_reads == 0:
                    continue

                print "Found %d region" %(relevant_region)
                relevant_region += 1

                # Compute the insert length of each read
                for read_pair_id, read_pair in paired_reads.iteritems():
                    if len(read_pair) != 2:
                        # Skip non-paired reads
                        continue
                    
                    left_read, right_read = read_pair
                    insert_len = right_read.pos - left_read.pos + 1

                    if insert_len > 0:
                        insert_lengths.append(insert_len)
                    else:
                        print "Negative or zero insert length ignored..."

    # Output results to file
    output_file = open(output_filename, 'w')
    insert_length_str = "\n".join(map(str, insert_lengths))
    output_file.write(insert_length_str)
    output_file.close()
                    
    t2 = time.time()
    print "Insert length computation took %.2f seconds." %(t2 - t1)
    

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-insert-len", dest="compute_insert_len", nargs=3, default=None,
                      help="Compute insert length for given sample. Takes as input "
                      "a sorted, indexed BAM file with headers, a GFF file with gene models, "
                      "and an output directory.")
    parser.add_option("--min-exon-size", dest="min_exon_size", nargs=1, type="int", default=500,
                      help="Minimum size of constitutive exon (in nucleotides) that should be used "
                      "in the computation.  Default is 500 bp.")
    (options, args) = parser.parse_args()

    if options.compute_insert_len != None:
        bam_filename = os.path.abspath(os.path.expanduser(options.compute_insert_len[0]))
        gff_filename = os.path.abspath(os.path.expanduser(options.compute_insert_len[1]))
        output_dir = os.path.abspath(os.path.expanduser(options.compute_insert_len[2]))

        compute_insert_len(bam_filename, gff_filename, output_dir, options.min_exon_size)


if __name__ == "__main__":
    main()
