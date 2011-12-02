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
import exon_utils


def load_insert_dist(insert_dist_filename):
    """
    Read insert length distribution.
    """
    insert_dist_file = open(insert_dist_filename, "r")
    insert_dist = array([int(line.strip()) \
                         for line in insert_dist_file])
    return insert_dist


def bedtools_map_bam_to_bed(bam_filename, gff_intervals_filename):
    """
    Map BAM file to GFF intervals and return the result as a
    BED file.

    Returns a stream to a BED file with the results
    """
    bedtools_cmd = "intersectBed -abam %s -b %s -wa -wb -bed -f 1" \
                   %(bam_filename, gff_intervals_filename)

    print "Executing: %s" %(bedtools_cmd)

    if (not os.path.isfile(bam_filename)) or \
       (not os.path.isfile(gff_intervals_filename)):
        raise Exception, "Error: %s or %s do not exist." %(bam_filename,
                                                           gff_intervals_filename)
    bed_stream = os.popen(bedtools_cmd)
    return bed_stream 

    
def compute_insert_len(bams_to_process, gff_filename, output_dir,
                       min_exon_size):
    """
    Compute insert length distribution and output it to the given
    directory.

    Arguments:

    - bams_to_process: a list of BAM files to process
    - gff_filename: GFF with gene models
    """
    bams_str = "\n  ".join(bams_to_process)
    num_bams = len(bams_to_process)
    print "Computing insert length distribution of %d files:\n  %s" \
          %(num_bams, bams_str)
    print "  - Using gene models from: %s" %(gff_filename)
    print "  - Outputting to: %s" %(output_dir)
    print "  - Minimum exon size used: %d" %(min_exon_size)

    if not os.path.isdir(output_dir):
        print "Making directory: %s" %(output_dir)
        os.makedirs(output_dir)

    # Get the constitutive exons that meet the size requirement
    const_exons, gff_exons_filename = exon_utils.get_const_exons_by_gene(gff_filename, output_dir,
                                                                         min_size=min_exon_size)

    for bam_filename in bams_to_process:
        output_filename = os.path.join(output_dir,
                                       "%s.insert_len" %(os.path.basename(bam_filename)))

        # Load BAM file with reads
        #bamfile = sam_utils.load_bam_reads(bam_filename)

        t1 = time.time()
        bed_mapped_reads = bedtools_map_bam_to_bed(bam_filename, gff_exons_filename)
        t2 = time.time()

        print "Fetching of reads in constitutive exons took: %.2f seconds" \
              %(t1 - t2)

        for mapped_read in bed_mapped_reads:
            print "=> ", mapped_read

        insert_lengths = []

        t1 = time.time()

        relevant_region = 0
    
        # paired_reads = sam_utils.pair_sam_reads(exon_reads)
        # num_paired_reads = len(paired_reads)

        # if num_paired_reads == 0:
        #     continue
        # for read_pair_id, read_pair in paired_reads.iteritems():
        #     if len(read_pair) != 2:
        #         # Skip non-paired reads
        #         continue

        #     left_read, right_read = read_pair
        #     insert_len = right_read.pos - left_read.pos + 1

        #     if insert_len > 0:
        #         insert_lengths.append(insert_len)
        #     else:
        #         print "Negative or zero insert length ignored..."

        # Output results to file
        #output_file = open(output_filename, 'w')
        #insert_length_str = "\n".join(map(str, insert_lengths))
        #output_file.write(insert_length_str)
        #output_file.close()

        t2 = time.time()
        print "Insert length computation took %.2f seconds." %(t2 - t1)


# def pair_reads_from_bed_intervals(bed_stream):
#     """
#     Match up read mates with each other, indexed by the BED interval
#     that they fall in.

#     Return a dictionary of BED region mapping to a set of read pairs.

#     Arguments:

#     - bed_filename: file with BED reads and the region they map to.

#     Returns. 
#     """
#     return

# def compute_insert_len(bam_filename, gff_filename, output_dir,
#                        min_exon_size):
#     """
#     Compute insert length distribution and output it to the given
#     directory.
#     """
#     print "Computing insert length distribution of %s" %(bam_filename)
#     print "  - Using gene models from: %s" %(gff_filename)
#     print "  - Outputting to: %s" %(output_dir)
#     print "  - Minimum exon size used: %d" %(min_exon_size)

#     if not os.path.isdir(output_dir):
#         print "Making directory: %s" %(output_dir)
#         os.makedirs(output_dir)

#     output_filename = os.path.join(output_dir,
#                                    "%s.insert_len" %(os.path.basename(bam_filename)))

#     # Load BAM file with reads
#     bamfile = sam_utils.load_bam_reads(bam_filename)
    
#     # Load the genes from the GFF
#     print "Loading genes from GFF..."
#     t1 = time.time()
#     gff_genes = gene_utils.load_genes_from_gff(gff_filename)
#     t2 = time.time()
#     print "  - Loading genes from GFF took %.2f seconds" %(t2 - t1)

#     insert_lengths = []

#     t1 = time.time()

#     relevant_region = 0
    
#     for gene_id, gene_info in gff_genes.iteritems():
#         gene_obj = gene_info["gene_object"]

#         # Get all the constitutive parts
#         const_parts = gene_obj.get_const_parts()

#         chrom = gene_obj.chrom

#         # Consider only the large constitutive parts
#         for part in const_parts:
#             if part.len >= min_exon_size:
#                 # Get all the reads that land in the coordinates of the exon
#                 try:
#                     exon_reads = bamfile.fetch(chrom, part.start, part.end)
#                 except ValueError:
#                     print "Could not fetch from region: ", chrom, part.start, part.end
#                     continue

#                 # Pair all the paired-end reads that land there
#                 paired_reads = sam_utils.pair_sam_reads(exon_reads)
#                 num_paired_reads = len(paired_reads)

#                 if num_paired_reads == 0:
#                     continue

#                 print "Found %d region" %(relevant_region)
#                 relevant_region += 1

#                 # Compute the insert length of each read
#                 for read_pair_id, read_pair in paired_reads.iteritems():
#                     if len(read_pair) != 2:
#                         # Skip non-paired reads
#                         continue
                    
#                     left_read, right_read = read_pair
#                     insert_len = right_read.pos - left_read.pos + 1

#                     if insert_len > 0:
#                         insert_lengths.append(insert_len)
#                     else:
#                         print "Negative or zero insert length ignored..."

#     # Output results to file
#     output_file = open(output_filename, 'w')
#     insert_length_str = "\n".join(map(str, insert_lengths))
#     output_file.write(insert_length_str)
#     output_file.close()
                    
#     t2 = time.time()
#     print "Insert length computation took %.2f seconds." %(t2 - t1)
    

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-insert-len", dest="compute_insert_len", nargs=3, default=None,
                      help="Compute insert length for given sample. Takes as input "
                      "(1) a comma-separated list of sorted, indexed BAM files with headers "
                      "(or a single BAM filename), (2) a GFF file with gene models, "
                      "and (3) an output directory.")
    parser.add_option("--min-exon-size", dest="min_exon_size", nargs=1, type="int", default=500,
                      help="Minimum size of constitutive exon (in nucleotides) that should be used "
                      "in the computation. Default is 500 bp.")
    (options, args) = parser.parse_args()

    if options.compute_insert_len != None:
        bams_to_process = [os.path.abspath(os.path.expanduser(f)) for f in \
                           options.compute_insert_len[0].split(",")]
        gff_filename = os.path.abspath(os.path.expanduser(options.compute_insert_len[1]))
        output_dir = os.path.abspath(os.path.expanduser(options.compute_insert_len[2]))

        compute_insert_len(bams_to_process, gff_filename, output_dir,
                           options.min_exon_size)


if __name__ == "__main__":
    main()
