#!/usr/bin/env python
##
## Script to build an indexed representation of a GFF file for efficient
## retrieval of genes
##
import os
import time
import glob
import shelve

from collections import defaultdict

import misopy
import misopy.gff_utils as gff_utils
import misopy.pickle_utils as pickle_utils
import misopy.Gene as gene_utils

def serialize_genes(gff_genes, output_dir):
    """
    Output genes into pickle files by chromosome, by gene.
    """
    genes_by_chrom = defaultdict(dict)

    # Split up genes by chromosome 
    for gene_id, gene_info in gff_genes.iteritems():
        gene_obj = gene_info["gene_object"]
        gene_hierarchy = gene_info["hierarchy"]
        genes_by_chrom[gene_obj.chrom][gene_id] = {'gene_object': gene_obj,
                                                   'hierarchy': gene_hierarchy}

    # Mapping from gene IDs to pickled filename
    gene_id_to_filename = {}
                                                   
    # Serialize all the genes in each chromosome into their
    # own directory
    for chrom, chrom_genes in genes_by_chrom.iteritems():
        if chrom.startswith("chr"):
            chrom_dir_name = chrom
        else:
            chrom_dir_name = "chr%s" %(str(chrom))

        # Make directory for chromosome if it doesn't already exist
        chrom_dir = os.path.join(output_dir, chrom_dir_name)
        if not os.path.isdir(chrom_dir):
            print "Making directory: %s" %(chrom_dir)
            os.makedirs(chrom_dir)

        t1 = time.time()
        # Serialize each gene into a separate file
        num_genes = len(genes_by_chrom[chrom])
        
        for gene_id, gene_info in genes_by_chrom[chrom].iteritems():
            gene_filename = os.path.abspath(os.path.join(chrom_dir,
                                                         "%s.pickle" %(gene_id)))
            pickle_utils.write_pickled_file({gene_id: genes_by_chrom[chrom][gene_id]},
                                            gene_filename)
            # Record what filename was associated with this gene ID
            gene_id_to_filename[gene_id] = gene_filename
            
        t2 = time.time()
        print "  - Chromosome serialization took %.2f seconds" %(t2 - t1)

    # Shelve the mapping from gene ids to filenames
    shelved_filename = os.path.join(output_dir, "genes_to_filenames.shelve")
    shelved_data = shelve.open(shelved_filename)
    
    for k, v in gene_id_to_filename.iteritems():
        shelved_data[k] = v
    shelved_data.close()
    
        
def index_gff(gff_filename, output_dir):
    """
    Index the given GFF and placed the indexed representation
    in the output directory.
    """
    print "Indexing GFF..."

    # First check that the GFF is not already indexed
    indexed_files = glob.glob(os.path.join(output_dir, "chr*"))
    if len(indexed_files) >= 1:
        print "%s appears to already be indexed. Aborting." %(gff_filename)
        return
    
    print "  - GFF: %s" %(gff_filename)
    print "  - Outputting to: %s" %(output_dir)
    overall_t1 = time.time()
    t1 = time.time()
    gff_genes = gene_utils.load_genes_from_gff(gff_filename)
    t2 = time.time()
    print "  - Loading of genes from GFF took %.2f seconds" %(t2 - t1)

    t1 = time.time()
#    pickle_filename = os.path.join(output_dir,
#                                   "%s.pickle" %(os.path.basename(os.path.abspath(gff_filename))))
#    pickle_utils.write_pickled_file(gff_genes, pickle_filename)
    serialize_genes(gff_genes, output_dir)
    t2 = time.time()
    print "  - Serialization of genes from GFF took %.2f seconds" %(t2 - t1)
    overall_t2 = time.time()
    print "Indexing of GFF took %.2f seconds." %(overall_t2 - overall_t1)

    
def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--index", dest="index_gff", nargs=2, default=None,
                      help="Index the given GFF. Takes as arguments as GFF filename "
                      "and an output directory.")
    (options, args) = parser.parse_args()

    if options.index_gff != None:
        gff_filename = os.path.abspath(os.path.expanduser(options.index_gff[0]))
        output_dir = os.path.abspath(os.path.expanduser(options.index_gff[1]))

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        index_gff(gff_filename, output_dir)


if __name__ == '__main__':
    main()
