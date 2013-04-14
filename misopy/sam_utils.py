##
## Utilities for handling SAM/BAM reads
##

from collections import defaultdict

import misopy
from misopy.Gene import load_genes_from_gff

import os
import time
import pysam
import binascii
import ctypes

import numpy as np
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


def single_end_read_to_isoforms(read, gene, read_len, overhang_len=1):
    """
    Align single-end SAM read to gene's isoforms.
    """
    start = read.pos + 1
    end = cigar_to_end_coord(start, read.cigar)    
    
    assert(start < end)

    alignment, isoform_coords = \
        gene.align_read_to_isoforms_with_cigar(read.cigar, start, end, read_len,
                                               overhang_len)
    return alignment


def paired_read_to_isoforms(paired_read, gene, read_len,
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
#    print "LEFT read start,end: ", left_start, left_end
#    print "RIGHT read start,end: ", right_start, right_end
    if left_start > right_start:
        return None, None
    else:
        alignment, frag_lens = gene.align_read_pair_with_cigar(
            left_read.cigar, left_start, left_end, 
            right_read.cigar, right_start, right_end,
            read_len=read_len, overhang=overhang_len)
        
    
#    assert(left_start < right_start), "Reads scrambled?"
#    assert(left_end < right_start), "Reads scrambled: left=(%d, %d), right=(%d, %d)" \
#                    %(left_start, left_end, right_start, right_end)
    
#    print "Read: ", (left_start, left_end), " - ", (right_start, right_end)
#    print "  - Alignment: ", alignment, " frag_lens: ", frag_lens
#    print "  - Sequences: "
#    print "  - %s\t%s" %(left_read.seq, right_read.seq)
    
    return alignment, frag_lens

# def paired_read_to_isoforms(paired_read, gene, read_len=36,
#                             overhang_len=1):
#     """
#     Align paired-end SAM read to gene's isoforms.
#     """
#     left_read = paired_read[0]
#     right_read = paired_read[1]

#     # Convert to 1-based coordinates
#     left_start = left_read.pos + 1
#     right_start = right_read.pos + 1

#     # Get end coordinates of each read
#     left_end = cigar_to_end_coord(left_start, left_read.cigar)

#     assert(left_start < left_end)
    
#     right_end = cigar_to_end_coord(right_start, right_read.cigar)

#     assert(right_start < right_end)

#     alignment = None
#     frag_lens = None
    
#     # Throw out reads that posit a zero or less than zero fragment
#     # length
#     if left_start > right_start or left_end > right_start:
#         return None, None
#     else:
#         alignment, frag_lens = gene.align_read_pair(left_start, left_end,
#                                                     right_start, right_end,
#                                                     read_len=read_len,
#                                                     overhang=overhang_len)
        
    
# #    assert(left_start < right_start), "Reads scrambled?"
# #    assert(left_end < right_start), "Reads scrambled: left=(%d, %d), right=(%d, %d)" \
# #                    %(left_start, left_end, right_start, right_end)
    
# #    print "Read: ", (left_start, left_end), " - ", (right_start, right_end)
# #    print "  - Alignment: ", alignment, " frag_lens: ", frag_lens
# #    print "  - Sequences: "
# #    print "  - %s\t%s" %(left_read.seq, right_read.seq)
    
#     return alignment, frag_lens
                         
def load_bam_reads(bam_filename,
                   template=None):
    """
    Load a set of indexed BAM reads.
    """
    print "Loading BAM filename from: %s" %(bam_filename)
    bam_filename = os.path.abspath(os.path.expanduser(bam_filename))
    bamfile = pysam.Samfile(bam_filename, "rb",
                            template=template)
    return bamfile


def fetch_bam_reads_in_gene(bamfile, chrom, start, end, gene=None):
    """
    Align BAM reads to the gene model.
    """
    gene_reads = []

    if chrom in bamfile.references:
        pass
    else:
        chrom_parts = chrom.split("chr")
        if len(chrom_parts) <= 1:
            chrom = chrom_parts[0]
        else:
            chrom = chrom_parts[1]

    try:
        gene_reads = bamfile.fetch(chrom, start, end)
    except ValueError:
        print "Cannot fetch reads in region: %s:%d-%d" %(chrom,
                                                         start,
                                                         end)
    except AssertionError:
        print "AssertionError in region: %s:%d-%d" %(chrom,
                                                     start,
                                                     end)
        print "  - Check that your BAM file is indexed!"
    return gene_reads


def flag_to_strand(flag):
    """
    Takes integer flag as argument.
    Returns strand ('+' or '-') from flag.
    """
    if flag == 0 or not (int(bin(flag)[-5]) & 1):
        return "+"
    return "-"


def strip_mate_id(read_name):
    """
    Strip canonical mate IDs for paired end reads, e.g.
    
    #1, #2

    or:

    /1, /2
    """
    if read_name.endswith("/1") or read_name.endswith("/2") or \
       read_name.endswith("#1") or read_name.endswith("#2"):
        read_name = read_name[0:-3]
    return read_name

    
def pair_sam_reads(samfile, filter_reads=True,
                   return_unpaired=False):
    """
    Pair reads from a SAM file together.
    """
    paired_reads = defaultdict(list)
    unpaired_reads = {}

    for read in samfile:
        curr_name = read.qname

        # Strip canonical mate IDs 
        curr_name = strip_mate_id(curr_name)
        
        if filter_reads:
            # Skip reads that failed QC or are unmapped
            if read.is_qcfail or read.is_unmapped or \
               read.mate_is_unmapped or (not read.is_paired):
                unpaired_reads[curr_name] = read
                continue
        paired_reads[curr_name].append(read)

    to_delete = []

    num_unpaired = 0
    num_total = 0

    for read_name, read in paired_reads.iteritems():
        if len(read) != 2:
            unpaired_reads[read_name] = read
            num_unpaired += 1
            # Delete unpaired reads
            to_delete.append(read_name)
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
            print "WARNING: %s left mate starts later than right "\
                  "mate" %(left_read.qname)
        num_total += 1

    # Delete reads that are on the same strand
    for del_key in to_delete:
        del paired_reads[del_key]

    print "Filtered out %d read pairs that were on same strand." \
        %(len(to_delete))
    print "Filtered out %d reads that had no paired mate." \
        %(num_unpaired)
    print "  - Total read pairs: %d" %(num_total)

    if not return_unpaired:
        return paired_reads
    else:
        return paired_reads, unpaired_reads


# Global variable containing CIGAR types for conversion
CIGAR_TYPES = ('M', 'I', 'D', 'N', 'S', 'H', 'P')

def sam_cigar_to_str(sam_cigar):
    """
    Convert pysam CIGAR list to string format.
    """
    # First element in sam CIGAR list is the CIGAR type
    # (e.g. match or insertion) and the second is
    # the number of nucleotides
    #cigar_str = "".join(["%d%s" %(c[1], CIGAR_TYPES[c[0]]) \
    #                     for c in sam_cigar])
    #### OPTIMIZED VERSION
    cigar_str = ""
    for c in sam_cigar:
        cigar_str += "%d%s" %(c[1], CIGAR_TYPES[c[0]])
    return cigar_str


def read_matches_strand(read,
                        target_strand,
                        strand_rule,
                        paired_end=None):
    """
    Check if a read matches strand.

    - target_strand: the annotation strand ('+' or '-')
    - strand_rule: the strand rule, i.e.
      ('fr-unstranded', 'fr-firststrand', or 'fr-secondstrand')
    """
    if strand_rule == "fr-unstranded":
        return True
    matches = False
    if paired_end is not None:
        # Paired-end reads
        read1, read2 = read
        if strand_rule == "fr-firststrand":
            # fr-firststrand: means that the *second* of the mates
            # must match the strand
            matches = (flag_to_strand(read2.flag) == target_strand)
        elif strand_rule == "fr-secondstrand":
            # fr-secondstrand: means that the *first* of the mates
            # must match the strand
            matches = (flag_to_strand(read1.flag) == target_strand)
        else:
            raise Exception, "Unknown strandedness rule."
    else:
        # Single-end reads
        if strand_rule == "fr-firststrand":
            # fr-firststrand: We sequence the first read only, so it must
            # *NOT* match the target strand
            matches = (flag_to_strand(read.flag) != target_strand)
        elif strand_rule == "fr-secondstrand":
            # fr-secondstrand: We only sequence the first read, which 
            # is supposed to match the target strand
            matches = (flag_to_strand(read.flag) == target_strand)
        else:
            raise Exception, "Unknown strandedness rule."
    return matches 


def sam_parse_reads(samfile,
                    paired_end=False,
                    strand_rule=None,
                    target_strand=None):
    """
    Parse the SAM reads. If paired-end, pair up the mates
    together.

    Also forces the strandedness convention, discarding
    reads that do not match the correct strand.

    - strand_rule: specifies the strandedness convention. Can be
      'fr-unstranded', 'fr-firststrand' or 'fr-secondstrand'.
    - target_strand: specifies the strand to match, i.e. the
      annotation strand. Can be '+' or '-'.
    """
    read_positions = []
    read_cigars = []
    num_reads = 0

    check_strand = True
    # Determine if we need to check strandedness of reads.
    # If we're given an unstranded convention, or if we're
    # not given a target strand, then assume that there's
    # no need to check strandedness.
    if (strand_rule is None) or \
       (strand_rule is "fr-unstranded") or \
       (target_strand is None):
        # No need to check strand
        check_strand = False

    # Track number of reads discarded due to strand
    # violations, if strand-specific
    num_strand_discarded = 0
    if paired_end:
        # Pair up the reads 
        paired_reads = pair_sam_reads(samfile)

        # Process reads into format required by fastmiso
        # MISO C engine requires pairs to follow each other in order.
        # Unpaired reads are not supported.
        for read_id, read_info in paired_reads.iteritems():
            if check_strand:
                # Check strand
                if not read_matches_strand(read_info,
                                           target_strand,
                                           strand_rule,
                                           paired_end=paired_end):
                    # Skip reads that don't match strand
                    num_strand_discarded += 1
                    continue
            read1, read2 = read_info
            # Read positions and cigar strings are collected
            read_positions.append(int(read1.pos))
            read_positions.append(int(read2.pos))
            read_cigars.append(sam_cigar_to_str(read1.cigar))
            read_cigars.append(sam_cigar_to_str(read2.cigar))
            num_reads += 1
    else:
        # Single-end
        for read in samfile:
            if check_strand:
                if not read_matches_strand(read,
                                           target_strand,
                                           strand_rule,
                                           paired_end=paired_end):
                    # Skip reads that don't match strand
                    num_strand_discarded += 1
                    continue
            read_positions.append(int(read.pos))
            read_cigars.append(sam_cigar_to_str(read.cigar))
            num_reads += 1

    if check_strand:
        print "No. reads discarded due to strand violation: %d" \
            %(num_strand_discarded)

    reads = (tuple(read_positions),
             tuple(read_cigars))

    return reads, num_reads
    

def sam_pe_reads_to_isoforms(samfile, gene, read_len, overhang_len):
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
        
        alignment, frag_lens = paired_read_to_isoforms(read_pair, gene,
                                                       read_len, overhang_len)
        
        # Skip reads that are not consistent with any isoform
        if any(array(alignment) == 1):
            pe_reads.append([alignment, frag_lens])
            num_read_pairs += 1
        else:
#            print "read %s inconsistent with all isoforms" %(read_id)
            k += 1

    print "Filtered out %d reads that were not consistent with any isoform" %(k)
    return pe_reads, num_read_pairs


# def sam_se_reads_to_isoforms(samfile, gene, read_len,
#                              overhang_len):
#     """
#     Align single-end reads to gene.
#     """
#     num_reads = 0

#     alignments = []

#     num_skipped = 0
    
#     for read in samfile:
#         alignment = single_end_read_to_isoforms(read, gene, read_len,
#                                                 overhang_len)
#         if 1 in alignment:
#             # If the read aligns to at least one of the isoforms, keep it
#             alignments.append(alignment)
#             num_reads += 1
#         else:
#             num_skipped += 1

#     print "Skipped total of %d reads." %(num_skipped)
        
#     return alignments, num_reads

    

def sam_reads_to_isoforms(samfile, gene, read_len, overhang_len,
                          paired_end=None):
    """
    Align BAM reads to the gene model.
    """
    t1 = time.time()

    if paired_end != None:
        # Paired-end reads
        reads, num_reads = sam_pe_reads_to_isoforms(samfile, gene, read_len,
                                                    overhang_len)
    else:
        print "Single end"
        # Single-end reads
        #reads, num_reads = sam_se_reads_to_isoforms(samfile, gene, read_len,
        #                                            overhang_len)
        reads = se_reads_to_assignments(samfile, gene, read_len)
    
    t2 = time.time()
    print "%d reads" %(len(reads))
    print "Alignment to gene took %.2f seconds." %((t2 - t1))
    return reads


def get_isoform_cigars_at_pos(gene, genomic_start, read_len):
    """
    Compute cigar strings of isoforms. Given a genomic start position and
    a read length, compute the CIGAR strings that would arise
    from a read at that position for each isoform.

    Return a vector of CIGAR strings of length number of isoforms.
    """
    for isoform in gene.isoforms:
        # Convert the genomic start to isoform coordinate
        isoform_start = isoform.part_coord_to_isoform(genomic_start)
        # Compute CIGAR string for given read length


def read_to_isoform_assignment(gene, read, read_len):
    """
    Return assignment matrix for each of the gene's isoforms.
    """
    read_cigar_str = sam_cigar_to_str(read.cigar)
    # Get the CIGAR strings for each isoform
    isoform_cigars = \
        [isoform.get_cigar_from_start(read.pos, read_len) \
         for isoform in gene.isoforms]
    assignments = tuple([1 if iso_cigar == read_cigar_str else 0 \
                         for iso_cigar in isoform_cigars])
    return assignments


def se_reads_to_assignments(reads, gene, read_len, paired_end=False):
    """
    Calculate a read assignment matrix for reads to isoforms.
    """
    # For each read, align it to isoforms using its CIGAR string
    assignments = \
        [read_to_isoform_assignment(gene, read, read_len) \
         for read in reads]
    return assignments
        


# int splicing_assignment_matrix(const splicing_gff_t *gff, size_t gene,
# 			       int readLength, int overHang, 
# 			       splicing_matrix_t *matrix) {
#   size_t noiso;
#   splicing_vector_int_t exstart, exend, exidx;
#   size_t genestart=VECTOR(gff->start)[ VECTOR(gff->genes)[gene] ];
#   size_t geneend=VECTOR(gff->end)[ VECTOR(gff->genes)[gene] ];
#   size_t i, elen;
#   size_t p=0, lastp=geneend - genestart - readLength + 1;
#   splicing_vector_int_t cigar, mp, mppos;
#   splicing_vector_int_t isoseq, isomatch;
#   splicing_i_assignmat_data_t mpdata = { &mp, &mppos };

#   if (overHang > 1) { 
#     SPLICING_ERROR("Overhang is not implemented in assignment matrix yet.",
# 		   SPLICING_UNIMPLEMENTED);
#   }

#   SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));

#   SPLICING_CHECK(splicing_vector_int_init(&cigar, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &cigar);
#   SPLICING_CHECK(splicing_vector_int_init(&exstart, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &exstart);
#   SPLICING_CHECK(splicing_vector_int_init(&exend, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &exend);
#   SPLICING_CHECK(splicing_vector_int_init(&exidx, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &exidx);
#   SPLICING_CHECK(splicing_gff_exon_start_end(gff, &exstart, &exend, &exidx,
# 					     gene));
#   elen=splicing_vector_int_size(&exstart);
#   for (i=0; i<elen; i++) {
#     VECTOR(exstart)[i] -= genestart;
#     VECTOR(exend)[i] -= genestart;
#   }
  
#   SPLICING_CHECK(splicing_numeric_cigar(&exstart, &exend, &exidx, noiso, 0,
# 					&cigar));

#   splicing_vector_int_destroy(&exidx);
#   splicing_vector_int_destroy(&exend);
#   splicing_vector_int_destroy(&exstart);
#   SPLICING_FINALLY_CLEAN(3);
  
#   SPLICING_CHECK(splicing_vector_int_init(&mp, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &mp);
#   SPLICING_CHECK(splicing_vector_int_init(&mppos, noiso));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &mppos);
#   SPLICING_CHECK(splicing_vector_int_init(&isoseq, noiso));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &isoseq);
#   SPLICING_CHECK(splicing_vector_int_init(&isomatch, noiso));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &isomatch);

#   SPLICING_CHECK(splicing_matrix_resize(matrix, noiso, 0));

#   while (p <= lastp) {
#     size_t pos=0, nextp, tmppos, ipos, l;

#     /* Calculate the initial CIGAR strings for all isoforms */
#     splicing_vector_int_clear(&mp);
#     SPLICING_CHECK(splicing_vector_int_push_back(&mp, 0)); l=1;
#     for (i=0; i<noiso; i++) {
#       int rl=readLength;
#       pos++;
#       VECTOR(mppos)[i]=l+1;
#       if (VECTOR(cigar)[pos] <= 0) { 
# 	VECTOR(mppos)[i]=0;
# 	while (VECTOR(cigar)[pos] != 0) { pos++; }
#       } else {
# 	while (VECTOR(cigar)[pos] != 0) {
# 	  if (VECTOR(cigar)[pos] < 0) {
# 	    SPLICING_CHECK(splicing_vector_int_push_back(&mp, 
# 					 VECTOR(cigar)[pos])); l++;
# 	    pos++;
# 	  } else if (rl <= VECTOR(cigar)[pos]) {
# 	    SPLICING_CHECK(splicing_vector_int_push_back(&mp, rl)); l++;
# 	    while (VECTOR(cigar)[pos] != 0) { pos++; }
# 	    rl = 0;
# 	    break;
# 	  } else {
# 	    SPLICING_CHECK(splicing_vector_int_push_back(&mp, 
# 					 VECTOR(cigar)[pos])); l++;
# 	    rl=rl-VECTOR(cigar)[pos];
# 	    pos++;
# 	  }
# 	}
# 	if (rl != 0) { VECTOR(mppos)[i] = 0; }
#       }
#       SPLICING_CHECK(splicing_vector_int_push_back(&mp, 0)); l++;
#     } /* i<noiso */

#     /* Calculate the next position to check */
#     nextp = lastp + 1 - p;
#     tmppos = 0;
#     for (i=0; i<noiso; i++) {
#       tmppos++;
#       ipos=tmppos;
#       if (VECTOR(cigar)[tmppos]==0) { continue; }
#       if (abs(VECTOR(cigar)[tmppos]) < nextp) { 
# 	nextp = abs(VECTOR(cigar)[tmppos]);
#       }
#       if (VECTOR(cigar)[tmppos] > 0) {
# 	size_t pnextp=nextp;
# 	int j=ipos;
# 	int rl2=readLength;
# 	while (VECTOR(cigar)[j] != 0) {
# 	  if (VECTOR(cigar)[j] >= rl2) { 
# 	    pnextp=VECTOR(cigar)[j]-rl2+1; 
# 	    break;
# 	  } else if (VECTOR(cigar)[j] > 0) {
# 	    rl2 = rl2-VECTOR(cigar)[j];
# 	  }
# 	  j++;
# 	}
# 	if (pnextp < nextp) { nextp = pnextp; }
#       }
#       while (VECTOR(cigar)[tmppos] != 0) { tmppos++; }
#     }

#     /* Record the assignment classes. 
#        First we sort the cigar prefixes; then we calculate the
#        assignment class of the current position. Finally, we check
#        whether this class is already in the result and then update/add
#        it. */
#     for (i=0; i<noiso; i++) { VECTOR(isoseq)[i]=i; }
#     splicing_qsort_r(VECTOR(isoseq), noiso, sizeof(int), &mpdata, 
# 		     splicing_i_assignmat_cmp);
    
#     for (i=0; i<noiso && VECTOR(mppos)[ VECTOR(isoseq)[i] ]==0; i++) ;
#     while (i<noiso) {
#       int col, j;
#       i++;
#       splicing_vector_int_null(&isomatch);
#       VECTOR(isomatch)[ VECTOR(isoseq)[i-1] ] = 1;
#       while (i < noiso &&  ! splicing_i_assignmat_cmp(&mpdata, 
# 					      &(VECTOR(isoseq)[i-1]),
# 					      &(VECTOR(isoseq)[i]))) {
# 	VECTOR(isomatch)[ VECTOR(isoseq)[i] ] = 1;
# 	i++;
#       }
#       SPLICING_CHECK(splicing_matrix_add_cols(matrix, 1));
#       col=splicing_matrix_ncol(matrix)-1;
#       for (j=0; j<noiso; j++) { 
# 	MATRIX(*matrix, j, col) = VECTOR(isomatch)[j] * nextp;
#       }
#     }

#     /* Update the CIGAR strings */
#     tmppos=0, ipos=0;
#     for (i=0; i<noiso; i++) {
#       VECTOR(cigar)[ipos] = 0;
#       tmppos++; ipos++;
#       if (VECTOR(cigar)[tmppos] != 0) {
# 	VECTOR(cigar)[ipos] = 
# 	  (VECTOR(cigar)[tmppos] < 0 ? -1 : 1) * 
# 	  (abs(VECTOR(cigar)[tmppos])-nextp);
# 	tmppos++;
# 	if (VECTOR(cigar)[ipos] !=0 ) ipos++;
#       }
#       while (VECTOR(cigar)[tmppos] != 0) {
# 	VECTOR(cigar)[ipos]=VECTOR(cigar)[tmppos];
# 	tmppos++;
# 	ipos++;
#       }      
#     }
#     VECTOR(cigar)[ipos] = 0;
    
#     p = p + nextp;


#   } /* p <= lastp */
  
#   splicing_vector_int_destroy(&isomatch);
#   splicing_vector_int_destroy(&isoseq);
#   splicing_vector_int_destroy(&mppos);
#   splicing_vector_int_destroy(&mp);
#   splicing_vector_int_destroy(&cigar);
#   SPLICING_FINALLY_CLEAN(5);

#   splicing_i_assignmat_simplify(matrix);
  
#   return 0;
# }

# /* Calculate the CIGAR strings for the isoforms in 'subset', 
#    from 'sp1' first, and then from 'sp2'. Concatenate everything 
#    in mp, and index it in mppos. */

# int splicing_i_get_mp(const splicing_gff_converter_t *converter, 
# 		      const splicing_vector_int_t *subset,
# 		      splicing_vector_int_t *mp,
# 		      splicing_vector_int_t *mppos,
# 		      int sp1, int sp2, int readLength) {
  
#   int i, n=splicing_vector_int_size(subset); 
#   int l=0;
#   const splicing_vector_int_t *exstart=&converter->exstart;
#   const splicing_vector_int_t *exend=&converter->exend;
#   const splicing_vector_int_t *exidx=&converter->exidx;

#   splicing_vector_int_clear(mp); 
#   splicing_vector_int_clear(mppos);

#   SPLICING_CHECK(splicing_vector_int_push_back(mp, 0)); l++;
#   for (i=0; i<n; i++) {
#     int iso, pos, pos2, rl;

#     /* ------------------------------------------------------------ */
#     /* SP1 */

#     iso=VECTOR(*subset)[i];
#     pos=VECTOR(*exidx)[iso];
#     pos2=VECTOR(*exidx)[iso+1];
#     rl=readLength;

#     SPLICING_CHECK(splicing_vector_int_push_back(mppos, l));

#     /* Search for first exon. */
#     for (; pos < pos2 && VECTOR(*exend)[pos] < sp1; pos++) ;
#     if (pos==pos2) { 
#       printf("Looking for %i in\n", sp1);
#       splicing_vector_int_print(exstart);
#       splicing_vector_int_print(exend);
#       splicing_vector_int_print(exidx);
#       printf("between pos %i and %i\n", VECTOR(*exidx)[iso], pos2);
#       abort();
#       SPLICING_ERROR("Internal splicing error", SPLICING_EINTERNAL);
#     }
    
#     /* Record string */
#     while (pos < pos2 && rl>0) {
#       int len;
#       if (sp1 > VECTOR(*exstart)[pos]) { 
# 	len=VECTOR(*exend)[pos] - sp1 + 1;
#       } else {
# 	len=VECTOR(*exend)[pos] - VECTOR(*exstart)[pos];
#       }
      
#       if (len>=rl) { 
# 	SPLICING_CHECK(splicing_vector_int_push_back(mp, rl)); l++;
# 	rl=0;
# 	break;
#       } else { 
# 	int gap=VECTOR(*exstart)[pos+1] - VECTOR(*exend)[pos];
# 	SPLICING_CHECK(splicing_vector_int_push_back(mp, len)); l++; 
# 	SPLICING_CHECK(splicing_vector_int_push_back(mp, -gap)); l++;
# 	rl -= len;
# 	pos++;
#       }
#     }

#     /* ------------------------------------------------------------ */
#     /* SP2 */

#     pos=VECTOR(*exidx)[iso];
#     rl=readLength;

#     /* Search for first exon. */
#     for (; pos < pos2 && VECTOR(*exend)[pos] < sp2; pos++) ;
#     if (pos==pos2) { 
#       SPLICING_ERROR("Internal splicing error", SPLICING_EINTERNAL);
#     }
    
#     /* Record string */
#     while (pos < pos2 && rl>0) {
#       int len;
#       if (sp2 > VECTOR(*exstart)[pos]) { 
# 	len=VECTOR(*exend)[pos] - sp2 + 1;
#       } else {
# 	len=VECTOR(*exend)[pos] - VECTOR(*exstart)[pos];
#       }
      
#       if (len>=rl) { 
# 	SPLICING_CHECK(splicing_vector_int_push_back(mp, rl)); l++;
# 	rl=0;
# 	break;
#       } else { 
# 	int gap=VECTOR(*exstart)[pos+1] - VECTOR(*exend)[pos];
# 	SPLICING_CHECK(splicing_vector_int_push_back(mp, len)); l++; 
# 	SPLICING_CHECK(splicing_vector_int_push_back(mp, -gap)); l++;
# 	rl -= len;
# 	pos++;
#       }
#     }

#     SPLICING_CHECK(splicing_vector_int_push_back(mp, 0)); l++;
#   }
  
#   return 0;
# }






# int splicing_matchIso(const splicing_gff_t *gff, int gene, 
# 		      const splicing_vector_int_t *position, 
# 		      const char **cigarstr, int overHang, int readLength,
# 		      splicing_matrix_t *result) {

#   int noreads=splicing_vector_int_size(position);
#   int r, i;
#   splicing_vector_int_t exstart, exend, exidx, cigar, cigaridx, cigarlength;
#   size_t noiso;

#   if (overHang==0) { overHang=1; }
#   if (overHang < 1) {
#     SPLICING_ERROR("Overhang length invalid. Must be positive", 
# 		   SPLICING_EINVAL);
#   }
#   if (readLength < 0) { 
#     SPLICING_ERROR("Read length cannot be negative", SPLICING_EINVAL);
#   }

#   SPLICING_CHECK(splicing_gff_noiso_one(gff, gene, &noiso));
#   SPLICING_CHECK(splicing_vector_int_init(&exstart, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &exstart);
#   SPLICING_CHECK(splicing_vector_int_init(&exend, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &exend);
#   SPLICING_CHECK(splicing_vector_int_init(&exidx, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &exidx);
#   SPLICING_CHECK(splicing_gff_exon_start_end(gff, &exstart, &exend, 
# 					     &exidx, gene));

#   SPLICING_CHECK(splicing_vector_int_init(&cigar, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &cigar);
#   SPLICING_CHECK(splicing_vector_int_init(&cigaridx, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &cigaridx);
#   SPLICING_CHECK(splicing_vector_int_init(&cigarlength, 0));
#   SPLICING_FINALLY(splicing_vector_int_destroy, &cigarlength);
#   SPLICING_CHECK(splicing_parse_cigar(cigarstr, noreads, &cigar, &cigaridx, 
# 				      &cigarlength, readLength));

#   SPLICING_CHECK(splicing_matrix_resize(result, noiso, noreads));

#   for (r=0; r<noreads; r++) {
#     int *mycig=VECTOR(cigar) + VECTOR(cigaridx)[r];
#     int nocig=VECTOR(cigaridx)[r+1] - VECTOR(cigaridx)[r];
#     int len =VECTOR(cigarlength)[r];

#     /* If the read length is specified, then filter out the reads
#        that are shorter */
#     if (len < readLength) { 
#       for (i=0; i<noiso; i++) { MATRIX(*result, i, r)=0; }

#     /* We can also filter out the reads that do not satisfy the overhang 
#        constraint here. We assume that the CIGAR string starts and ends 
#        with a match (of non-zero length). */
#     } else if (mycig[0] < overHang || mycig[nocig-1] < overHang) {
#       for (i=0; i<noiso; i++) { MATRIX(*result, i, r)=0; }
#     } else {

#       for (i=0; i<noiso; i++) {
# 	int c, pos=VECTOR(*position)[r];
# 	int ex=VECTOR(exidx)[i];
	
# 	/* Look for the exon where the read starts */
# 	while (ex < VECTOR(exidx)[i+1] &&
# 	       (pos < VECTOR(exstart)[ex] || VECTOR(exend)[ex] < pos)) {
# 	  ex++;
# 	}
# 	if (ex >= VECTOR(exidx)[i+1]) { MATRIX(*result, i, r)=0; continue; }
	
# 	/* Got it, match cigar string to exons */
# 	MATRIX(*result, i, r)=1;
# 	for (c=0; c<nocig; c++) {
# 	  if (mycig[c] > 0) { /* exon */
# 	    if (pos + mycig[c] - 1 > VECTOR(exend)[ex]) { 
# 	      MATRIX(*result, i, r)=0; break;
# 	    }
# 	    pos += mycig[c];
# 	  } else {	  	   /* intron */
# 	    if (pos != VECTOR(exend)[ex]+1) {
# 	      MATRIX(*result, i, r)=0; break; 
# 	    }
# 	    pos -= mycig[c];
# 	    ex += 1;
# 	    if (ex >= VECTOR(exidx)[i+1] || pos != VECTOR(exstart)[ex]) {
# 	      MATRIX(*result, i, r)=0; break;
# 	    }
# 	  }
# 	}
#       } /* i < noiso */
#     }   /* r < noreads */
#   }	/* if overhang os OK */

#   splicing_vector_int_destroy(&cigarlength);
#   splicing_vector_int_destroy(&cigaridx);
#   splicing_vector_int_destroy(&cigar);
#   splicing_vector_int_destroy(&exidx);
#   splicing_vector_int_destroy(&exend);
#   splicing_vector_int_destroy(&exstart);
#   SPLICING_FINALLY_CLEAN(6);
  
#   return 0;
# }

# int splicing_getMatchVector(const splicing_gff_t *gff, int gene,
# 			    int no_reads, const splicing_vector_int_t *position,
# 			    const char **cigarstr, int overHang, int readLength,
# 			    const splicing_matrix_t *matchMatrix,
# 			    const splicing_matrix_t *assMatrix,
# 			    splicing_vector_t *match) {

#   int no_classes=splicing_matrix_ncol(assMatrix);
#   int noiso=splicing_matrix_nrow(assMatrix);
#   int r;

#   SPLICING_CHECK(splicing_vector_resize(match, no_classes));

#   for (r=0; r<no_reads; r++) {
#     size_t cl, found;
#     for (cl=0, found=0; !found && cl < no_classes; cl++) {
#       size_t i;
#       for (i=0, found=1; i<noiso && found; i++) {
# 	double m1=MATRIX(*matchMatrix, i, r);
# 	double m2=MATRIX(*assMatrix, i, cl);
# 	found = (m1 > 0 && m2 > 0) || (m1 == 0 && m2 == 0);
#       }
#     }
#     if (found) { VECTOR(*match)[cl-1] += 1; }
#   }

#   return 0;
# }


def main():
    pass


if __name__ == "__main__":
    main()
