##
## Utilities for handling SAM/BAM reads
##

from collections import defaultdict
from Gene import load_genes_from_gff

import time
import pysam
import binascii

from numpy import array
from scipy import *

# def read_to_isoforms(alignment, gene):
#     """
#     Align SAM read to gene's isoforms.
#     """
#     ##
#     ## SAM format is zero-based, but GFF are 1-based, so add 1
#     ## to read position
#     ##
#     pass

def cigar_to_end_coord(start, cigar):
    """
    Compute the end coordinate based on the CIGAR string.

    Assumes the coordinate is 1-based.
    """
    #print start, cigar
    # Parse cigar part
    #for cigar_part in cigar:
    #    cigar_type, seq_len = cigar_part
    #    offset += seq_len
    offset = sum([cigar_part[1] for cigar_part in cigar])
    end = start + offset - 1
    return end

def single_end_read_to_isoforms(read, gene, read_len=36, overhang_len=1):
    """
    Align single-end SAM read to gene's isoforms.
    """
    start = read.pos + 1
    end = cigar_to_end_coord(start, read.cigar)
    
    assert(start < end)
    
    alignment, isoform_coords = gene.align_read_to_isoforms(start, end,
                                                            overhang=overhang_len,
                                                            read_len=read_len)
    return alignment
    

def paired_read_to_isoforms(paired_read, gene, read_len=36,
                            overhang_len=1):
    """
    Align paired-end SAM read to gene's isoforms.
    """
    left_read = paired_read[0]
    right_read = paired_read[1]

    # Convert to 1-based coordinates
    left_start = left_read.pos + 1
    right_start = right_read.pos + 1

    # Get end coordinates of each read
    left_end = cigar_to_end_coord(left_start, left_read.cigar)

    assert(left_start < left_end)
    
    right_end = cigar_to_end_coord(right_start, right_read.cigar)

    assert(right_start < right_end)

    alignment = None
    frag_lens = None
    
    # Throw out reads that posit a zero or less than zero fragment
    # length
    if left_start > right_start or left_end > right_start:
        return None, None
    else:
        alignment, frag_lens = gene.align_read_pair(left_start, left_end,
                                                    right_start, right_end,
                                                    read_len=read_len,
                                                    overhang=overhang_len)
        
    
#    assert(left_start < right_start), "Reads scrambled?"
#    assert(left_end < right_start), "Reads scrambled: left=(%d, %d), right=(%d, %d)" \
#                    %(left_start, left_end, right_start, right_end)
    
#    print "Read: ", (left_start, left_end), " - ", (right_start, right_end)
#    print "  - Alignment: ", alignment, " frag_lens: ", frag_lens
#    print "  - Sequences: "
#    print "  - %s\t%s" %(left_read.seq, right_read.seq)
    
    return alignment, frag_lens
                         

def load_bam_reads(bam_filename,
                   template=None):
    """
    Load a set of indexed BAM reads.
    """
    print "Loading BAM filename from: %s" %(bam_filename)
    t1 = time.time()
    bamfile = pysam.Samfile(bam_filename, "rb",
                            template=template)
    t2 = time.time()
    print "Loading took %.2f seconds" %(t2 - t1)
    return bamfile


def load_sam_reads(sam_filename, template=None):
    """
    Load a set of SAM reads.
    """
    print "Loading SAM filename from: %s" %(sam_filename)
    t1 = time.time()
    print "Using header: ", template
    samfile = pysam.Samfile(sam_filename, "r",
                            template=template)
    t2 = time.time()
    print "Loading took %.2f seconds" %(t2 - t1)
    return samfile

def fetch_bam_reads_in_region(bamfile, chrom, start, end, gene=None):
    """
    Align BAM reads to the gene model.
    """
    gene_reads = []
    if chrom in bamfile.references:
        pass
    else:
        chrom = chrom.split("chr")[1]
    try:
        gene_reads = bamfile.fetch(chrom, start, end)
    except ValueError:
        print "Cannot fetch reads in region: %s:%d-%d" %(chrom,
                                                         start,
                                                         end)
    return gene_reads

    
def fetch_bam_reads_in_gene(bamfile, chrom, start, end, gene=None):
    """
    Align BAM reads to the gene model.
    """
    gene_reads = []
    if chrom in bamfile.references:
        pass
    else:
        chrom = chrom.split("chr")[1]
    try:
        gene_reads = bamfile.fetch(chrom, start, end)
    except ValueError:
        print "Cannot fetch reads in region: %s:%d-%d" %(chrom,
                                                         start,
                                                         end)
    return gene_reads


def flag_to_strand(flag):
    """
    Takes integer flag as argument.
    Returns strand ('+' or '-') from flag.
    """
    if flag == 0 or not (int(bin(flag)[-5]) & 1):
        return "+"
    return "-"


def pair_sam_reads(samfile, filter_reads=True,
                   return_unpaired=False):
    """
    Pair reads from a SAM file together.
    """
#    print "Pairing SAM reads..."
    paired_reads = defaultdict(list)
    unpaired_reads = {}

    for read in samfile:
        if filter_reads:
            # Skip reads that failed QC or are unmapped
            if read.is_qcfail or read.is_unmapped or \
               read.mate_is_unmapped or (not read.is_paired):
                unpaired_reads[read_name] = read
                continue
        paired_reads[read.qname].append(read)

    to_delete = []

    num_unpaired = 0

    for read_name, read in paired_reads.iteritems():
        if len(read) != 2:
            unpaired_reads[read_name] = read
            num_unpaired += 1
            continue
        left_read, right_read = read[0], read[1]

        # Check that read mates are on opposite strands
        left_strand = flag_to_strand(left_read.flag)
        right_strand = flag_to_strand(right_read.flag)

        if left_strand == right_strand:
            # Skip read pairs that are on the same strand
            to_delete.append(read_name)
            continue
        
        if left_read.pos > right_read.pos:
            raise Exception, (left_read.qname, left_read.pos, right_read.pos)

    # Delete reads that are on the same strand
    for del_key in to_delete:
        del paired_reads[del_key]

    print "Filtered out %d read pairs that were on same strand." %(len(to_delete))
    print "Filtered out %d reads that had no paired mate." %(num_unpaired)

    if not return_unpaired:
        return paired_reads
    else:
        return paired_reads, unpaired_reads


def sam_pe_reads_to_isoforms(samfile, gene):
    """
    Align read pairs (from paired-end data set) to gene.

    Returns alignment of paired-end reads (with insert lengths)
    to gene and number of read pairs aligned.
    """
    paired_reads = pair_sam_reads(samfile)

    num_read_pairs = 0

    pe_reads = []

    k = 0

    for read_id, read_pair in paired_reads.iteritems():
        if len(read_pair) != 2:
            # Skip reads with no pair
            continue
        
        alignment, frag_lens = paired_read_to_isoforms(read_pair, gene)

        # Skip reads that are not consistent with any isoform
        if any(array(alignment) == 1):
            pe_reads.append([alignment, frag_lens])
            num_read_pairs += 1
        else:
            k += 1

    print "Filtered out %d reads that were not consistent with any isoform" %(k)

    return pe_reads, num_read_pairs


def sam_se_reads_to_isoforms(samfile, gene):
    """
    Align single-end reads to gene.
    """
    num_reads = 0

    alignments = []
    
    for read in samfile:
        alignment = single_end_read_to_isoforms(read, gene)
        if 1 in alignment:
            # If the read aligns to at least one of the isoforms, keep it
            alignments.append(alignment)
            num_reads += 1
        
    return alignments, num_reads


def sam_reads_to_isoforms(samfile, gene, paired_end=False):
    """
    Align BAM reads to the gene model.
    """
    print "Aligning reads gene..."
    t1 = time.time()

    if paired_end != None:
        # Paired-end reads
        reads, num_reads = sam_pe_reads_to_isoforms(samfile, gene)
    else:
        # Single-end reads
        reads, num_reads = sam_se_reads_to_isoforms(samfile, gene)
    
    t2 = time.time()
    print "Alignment to gene took %.2f seconds (%d reads)." %((t2 - t1),
                                                              num_reads)
    return reads


    
def GRIA1_example(sam_filename):
    """
    Test GRIA1 example.
    """
    sam_filename = 'sam-tests/brain/GRIA1.mm9.stim30.sam'
    samfile = load_sam_reads(sam_filename)
#    chrom = 'chr1'
#    start = 1
#    end = 4500
    
    gene = None

    sam_reads_to_isoforms(samfile, gene)

    #align_sam_reads_to_gene(samfile, chrom, start, end, gene)
    

def Atp2b1_example(sam_filename, use_bam=False, template=None):
    """
    Test Atp2b1 example in C2C12.
    """
    if use_bam:
        print "Using BAM"
        samfile = load_bam_reads(sam_filename, template=template)
    else:
        print "Using SAM"
        samfile = load_sam_reads(sam_filename, template=template)
    gene = load_genes_from_gff('gff/mm9.Atp2b1.gff')

    for isoform in gene.isoforms:
        print "=" * 10
        print "isoform parts: "
        for p in isoform.parts:
            print " part: ", p
    sam_reads_to_isoforms(samfile, gene)


def load_sam_as_bam(header_filename=None):
    """
    Convert a SAM format into BAM, after sorting it and indexing it.
    """
    return
    

def main():
    #bam_filename = 'sam-tests/skeletal-muscle/accepted_hits.n250000.bam'
    #sam_filename = 'sam-tests/brain/GRIA1.mm9.sam'
    sam_filename = 'sam-tests/muscle-c2c12/c2c12_cugbp1_kd.control_0d.Atp2b1.sam'
    bam_filename = 'sam-tests/muscle-c2c12/c2c12_cugbp1_kd.control_0d.Atp2b1.sorted.bam'
    mm9_header = pysam.Samfile('sam-tests/muscle-c2c12/mm9.genome.fasta.fai', "r")
#    Atp2b1_example(bam_filename, template=mm9_header,
#                   use_bam=True)

#    bam_filename = 'sam-tests/brain/Grin1.mm9.control.sorted.bam'
    bamfile = load_bam_reads(bam_filename, template=mm9_header)
    for r in bamfile:
        print r
        break

if __name__ == "__main__":
    main()
