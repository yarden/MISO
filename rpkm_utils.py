##
## Utilities for computing gene expression in units of RPKM
## (reads per kilobase per milliom mapped reads.)
##
## NOTE: Requires bedtools to be available.
##

def compute_rpkm(bam_filename, exons_bed_filename, output_dir):
    """
    Compute RPKM from a set of constitutive exons, given in BED file format.

    Arguments:

    - filename for indexed, sorted BAM file with reads
    - filename with constitutive exons
    - output directory to output RPKM estimates to 
    """
    return


def main():
    return


if __name__ == '__main__':
    main()
