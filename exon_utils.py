##
## Utilities for fetching exons from GFF files
##
import time
import os

from collections import defaultdict
import GFF

def is_exon_in_mRNA(gff_in, exon, mRNA):
    """
    Check if the exon is in mRNA based on genomic coordinates.

    Arguments:

    - exon: gff record of exon
    - mRNA: gff record of mRNA
    """
    mRNA_id = mRNA.get_id()
    if mRNA_id not in gff_in.exons_by_mRNA:
        return False
    
    for curr_exon in gff_in.exons_by_mRNA[mRNA.get_id()]:
        # Exon in mRNA if it has the same coordinates and strand
        # as one of the mRNA's exons
        if curr_exon.start == exon.start and \
           curr_exon.end == exon.end and \
           curr_exon.strand == exon.strand:
            return True
    return False
    

def get_const_exons_from_mRNA(gff_in, mRNAs, min_size=0):
    const_exons = []
    # Get first mRNA's exons
    gene_id = mRNAs[0].get_parent()
    mRNA_id = mRNAs[0].get_id()
    if mRNA_id not in gff_in.exons_by_mRNA:
        # If mRNA has no exons, skip
        return const_exons

    exons = gff_in.exons_by_mRNA[mRNA_id]

    for exon in exons:
        const_exon = True
        # Skip exons that don't meet size requirement
        exon_len = exon.end - exon.start + 1
        if exon_len < min_size: continue
        
        # Exon flagged constitutive unless it is missing
        # in one of the mRNAs
        for mRNA in mRNAs[1:]:
            curr_mRNA_id = mRNA.get_id()
            if not is_exon_in_mRNA(gff_in, exon, mRNA):
                const_exon = False
                break

        # If exon is constitutive, add the parent gene information
        # as a field
        exon.attributes['GeneParent'] = [gene_id]
        # Record exon if it is constitutive
        if const_exon: const_exons.append(exon)

    return const_exons

def get_bam_reads_in_exons(exons_filename, bam_filename):
    """
    Get the BAM reads that land inside the boundaries of a set of exons.

    Arguments:

    - exons filename (in GFF format)

    - BAM filename

    Returns results in BED format.

    Relies on intersectBed being available.
    """
    #intersectBed -abam reads.bam -b exons.gff -wb -f 1 | coverageBed -abam stdin -b exons.gff
    cmd = "intersectBed -abam %s -b %s -wb -f 1 | coverageBed -abam stdin -b %s -bed" \
          %(bam_filename, exons_filename, exons_filename)
    print "Getting BAM reads in exons..."
#    a = subprocess.Popen("intersectBed -abam test-output/sam-output/c2c12.Atp2b1.sorted.bam -b test-gff/Atp2b1.mm9.gff.const_exons -wb -f 1 -bed".split(), stdout=subprocess.PIPE,stderr=subprocess.PIPE)    
    return bam_reads


def output_exons_to_file(recs, output_filename,
                         output_format='gff'):
    """
    Output exons to file.

    Arguments:

    - records in gff format
    - filename to output results to
    """
    print "Outputting exons to file: %s" %(output_filename)

    if output_format == "gff":
        # Write file in GFF format
        output_file = open(output_filename, 'w')
        gff_writer = GFF.Writer(output_file)
        recs.reverse()
        gff_writer.write_recs(recs)
        output_file.close()
    elif output_format == "bed":
        # Write file in BED format
        raise Exception, "BED format unsupported."
    

def get_const_exons_by_gene(gff_filename, output_dir,
                            min_size=0, output_format='gff'):
    """
    Get consitutive exons from GFF file.

    Arguments:
    - gff_filename: GFF input filename
    - output_dir: output directory

    Optional arguments:

    - min_size: minimum exon size
    - output_format: gff or BED
    """
    print "Getting constitutive exons..."
    print "  - Input GFF: %s" %(gff_filename)
    print "  - Output dir: %s" %(output_dir)
    print "  - Output format: %s" %(output_format)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if min_size > 0:
        print "  - Including only exons greater than %d-bp" \
              %(min_size)

    t1 = time.time()
    gff_in = GFF.GFFDatabase(from_filename=gff_filename)

    const_exons_by_gene = []

    num_exons = 0

    for gene, mRNAs in gff_in.mRNAs_by_gene.iteritems():
        # For each gene, look at all mRNAs and return constitutive exon
        curr_const_exons = get_const_exons_from_mRNA(gff_in, mRNAs,
                                                     min_size=min_size)
        const_exons_by_gene.extend(curr_const_exons)
        num_exons += len(curr_const_exons)

    t2 = time.time()

    print "Constitutive exon retrieval took %.2f seconds (%d exons)." \
          %((t2 - t1), num_exons)

    output_filename = os.path.join(output_dir,
                                   "%s.const_exons" %(os.path.basename(gff_filename)))

    output_exons_to_file(const_exons_by_gene, output_filename,
                         output_format=output_format)
    return const_exons_by_gene, output_filename


def main():
    test_gff = '/home/yarden/MISO/gff-events/mm9/genes/Atp2b1.mm9.gff'
    output_dir = 'test-gff'
    const_exons_by_gene, output_filename = get_const_exons_by_gene(test_gff, output_dir)
    print "const_exons_by_gene: ", const_exons_by_gene, output_filename

if __name__ == '__main__':
    main()
