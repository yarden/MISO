#!/usr/bin/env python
##
## Interface for running MISO locally or on cluster
##

import glob
import time
import os
import sys

import misopy
from misopy.settings import Settings
from misopy.settings import miso_path as miso_settings_path
import misopy.hypothesis_test as ht
import misopy.as_events as as_events
import misopy.cluster_utils as cluster_utils
import misopy.sam_utils as sam_utils
from misopy.parse_csv import *
from misopy.json_utils import *
import misopy.miso_sampler as miso
import misopy.Gene as gene_utils
import misopy.gff_utils as gff_utils
from misopy.samples_utils import *

import numpy as np
np.seterr(all='ignore')

miso_path = os.path.dirname(os.path.abspath(__file__))
    
def compute_two_iso_psi(events_filename, event_type, output_dir, read_len,
                        overhang_len, num_sampler_iters=5000):
    """
    Compute Psi values for the given set of events and dump them in
    the given output directory.
    """
    print "Running MISO on: %s" %(events_filename)
    print "  - Output directory: %s" %(output_dir)
    
    # Load two-isoform events
    miso_events = as_events.MISOEvents(2, event_type, from_file=events_filename)
    num_events = len(miso_events.events)

    # Run sampler on each event
    for event_name, event in miso_events.events.iteritems():
	run_two_iso_event(event.label, event_type, miso_events, output_dir,
                          read_len, overhang_len, num_sampler_iters)
    print "Submitted total of %d events." %(num_events)

    # Summarize results
    summary_dir = os.path.join(output_dir, 'summary')
    if not os.path.isdir(summary_dir):
	os.makedirs(summary_dir)
    summary_filename = os.path.join(summary_dir, 'events_summary.miso')
    summarize_sampler_results(output_dir, summary_filename)


def run_two_iso_event(event_name, event_type, miso_events, output_dir, read_len,
                      overhang_len, num_sampler_iters=5000):
    """
    Run MISO on a given event.
    """
    # Load two-isoform events
    #miso_events = as_events.MISOEvents(2, event_type, from_file=events_filename)
    
    # Parse only the current event to a gene
    events_to_genes = miso_events.loaded_events_to_genes(single_event_name=event_name,
                                                         read_len=read_len,
                                                         overhang_len=overhang_len)

    # Find event
    event = miso_events.get_event(event_name)

    # Get gene corresponding to event and run sampler on it
    gene = events_to_genes[event_name]
    
    if event_type == 'SE' or event_type == 'RI':
	ni = event.num_inc
	ne = event.num_exc
	nb = event.num_common
    elif event_type == 'TandemUTR':
	ni = event.num_ext
	ne = 0
	nb = event.num_core
    elif event_type == 'AFE':
        ni = event.num_proximal_body + event.num_proximal_jxns
        ne = event.num_distal_body + event.num_distal_jxns
        nb = 0
        print "Event %s has %d proximal, %d distal reads" %(event.label,
                                                            ni,
                                                            ne)
    else:
	raise Exception, "Unsupported event type: %s" %(event_type)
    samples, cred_interval = miso.run_sampler_on_event(gene, ni, ne, nb,
                                                       read_len, overhang_len,
                                                       num_sampler_iters, output_dir)
    

# def run_multi_iso_event(isoforms_filename, reads_filename, event_type, output_dir,
#                         read_len, overhang_len, num_sampler_iters=5000, burn_in=500,
#                         lag=10):
#     gene = gene_utils.load_multi_isoform_gene(isoforms_filename)
#     reads = gene_utils.load_multi_isoform_reads(reads_filename)
    
#     num_isoforms = len(gene.isoforms)
#     hyperparameters = ones(num_isoforms)
#     proposal_diag = 0.05
#     sigma = miso.set_diag(zeros([num_isoforms-1, num_isoforms-1]),
#                           proposal_diag)
#     params = {'read_len': read_len,
#               'overhang_len': overhang_len,
#               'uniform_proposal': False,
#               'sigma_proposal': sigma}
    
#     sampler = miso.MISOSampler(params)
    
#     if not os.path.isdir(output_dir):
#         os.makedirs(output_dir)
#     output_filename = os.path.join(output_dir, os.path.basename(isoforms_filename).split(".")[0])
#     sampler_results = sampler.run_sampler(num_sampler_iters, reads, gene, hyperparameters,
#                                           params, output_filename, burn_in=burn_in, lag=lag)
#     if sampler_results != None:
#         percent_acceptance, sampled_psi, total_log_scores, kept_log_scores = sampler_results
#         cred_intervals = ht.compute_multi_iso_credible_intervals(sampled_psi)
#         return sampled_psi, cred_intervals
    
#     return None
    
    
def run_two_iso_on_cluster(miso_path, events_filename, event_type, psi_outdir,
                           read_len, overhang_len, chunk_jobs=False):
    """
    Run two-isoform MISO on cluster.

    - chunk_jobs: the number of jobs in each batch.  All jobs in a batch will be assigned to the same processor on
      the cluster.  When chunk_jobs is not specified, each event gets sent as a separate job.
    """
    print "Running two isoform MISO on cluster..."
    # Load two-isoform events
    miso_events = as_events.MISOEvents(2, event_type, from_file=events_filename)
    num_total_events = len(miso_events.events)
    delay_constant = 0.9

    if not chunk_jobs:
	event_batches = [miso_events.events]
    else:
	# Make sure we're chunking into more than one batch of jobs
	assert(chunk_jobs > 1)

        # Compute number of chunks we'd need to split all events to in order to get
	# 'chunk_jobs'-many events in a job
	chunk_jobs = int(round(num_total_events / float(chunk_jobs)))
	print "Splitting events into %d chunks..." %(chunk_jobs)
	event_names = miso_events.events.keys()
	event_batches = cluster_utils.chunk_list(event_names, chunk_jobs)
	print "  - Total of %d event batches." %(len(event_batches))

    batch_lens = [len(batch) for batch in event_batches]
    max_events_per_batch = max(batch_lens)
    queue_thresh = 50
    num_batches = len(event_batches)
    long_batch = 100
    
    if max_events_per_batch >= queue_thresh and max_events_per_batch <= long_batch:
	print "Longest batch contains more than %d jobs -- changing queue type to short" \
              %(queue_thresh)
	queue_type = 'short'
    else:
        print "Longest batch contains more than %d jobs -- changing queue type to long" \
              %(long_batch)
        queue_type = 'long'

    for event_batch in event_batches:
        # Compile set of commands that will be run in the same job
	miso_event_cmd_list = []
	num_jobs_per_batch = len(event_batch)
	print "Processing a batch of size %d events" %(num_jobs_per_batch)
	for event_name in event_batch:
	    miso_event_cmd = 'python %s --run-two-iso-event \"%s\" %s %s --event-type %s --read-len %d --overhang-len %d' \
			     %(os.path.join(miso_path, 'run_miso.py'),
			       event_name,
			       events_filename,
			       psi_outdir,
			       event_type,
			       read_len,
			       overhang_len)
	    miso_event_cmd_list.append(miso_event_cmd)
	# Execute events in batch
	miso_event_batch_cmd = "; ".join(miso_event_cmd_list)
	#print "Executing batched command list: ", miso_event_batch_cmd
	if num_batches > 1:
	    event_name += "_batch"
	cluster_utils.run_on_cluster(miso_event_batch_cmd, event_name, psi_outdir,
                                     queue_type=queue_type)
	# Add pause to allow cluster to process jobs
	time.sleep(delay_constant)
    # Parse all events into genes
    events_to_genes = miso_events.loaded_events_to_genes(read_len=read_len,
                                                         overhang_len=overhang_len)

def get_current_args():
    """
    Return the current arguments as a string.
    """
    return " ".join(sys.argv)

def get_curr_script_cmd():
    """
    Get the invocation of the current script (with its command line arguments) as a
    full command for use in a script.
    """
    return 'python ' + get_current_args()

def strip_option(cmd, option):
    """
    Strip given option from the given command line argument.
    """
    return "".join(cmd.split(option))

def get_bayes_factor_filenames(comparison_dir):
    """
    Given a comparison directory, return a list of all the filenames
    that end in .miso_bf and are in a _vs_ comparison sub-directory.
    """
    bf_files_matcher = os.path.join(comparison_dir,
                                    '*_vs_*',
                                    'bayes-factors',
                                    '*.miso_bf')
    comparison_filenames = glob.glob(bf_files_matcher)
    return comparison_filenames

def get_psi_info_by_sample(event_comparison_data, sample1_or_sample2):
    psi_info = {}
    
    for old_key in event_comparison_data.keys():
        if old_key.startswith(sample1_or_sample2):
            new_key = old_key.split("%s_" %(sample1_or_sample2))[1]
            psi_info[new_key] = event_comparison_data[old_key]

    return psi_info


##
## Multi-isoform interface
##
def compute_gene_psi(gene_ids, gff_index_filename, bam_filename, output_dir,
                     read_len, overhang_len, paired_end=None, event_type=None):
    """
    Run Psi at the Gene-level (for multi-isoform inference.)

    Arguments:

    - Set of gene IDs corresponding to gene IDs from the GFF
    - Indexed GFF filename describing the genes
    - BAM filename with the reads (must be sorted and indexed)
    - Output directory
    - Optional: Run in paired-end mode. Gives mean and standard deviation
      of fragment length distribution.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(gff_index_filename):
        print "Error: no such GFF file as %s" %(gff_index_filename)
        return

    num_genes = len(gene_ids)
    
    print "Computing Psi for %d genes..." %(num_genes)
    print "  - " + ", ".join(gene_ids)
    print "  - GFF filename: %s" %(gff_index_filename)
    print "  - BAM: %s" %(bam_filename)
    print "  - Outputting to: %s" %(output_dir)

    if paired_end:
        print "  - Paired-end mode: ", paired_end

    settings = Settings.get()
    settings_params = Settings.get_sampler_params()
    
    burn_in = settings_params["burn_in"]
    lag = settings_params["lag"]
    num_iters = settings_params["num_iters"]

    min_event_reads = Settings.get_min_event_reads()

    mean_frag_len = None
    frag_variance = None

    if paired_end:
        mean_frag_len = int(paired_end[0])
        frag_variance = power(int(paired_end[1]), 2)


    # Load the genes from the GFF
#    print "Loading genes from indexed GFF..."
#    t1 = time.time()
    gff_genes = gff_utils.load_indexed_gff_file(gff_index_filename)
#    t2 = time.time()
#    print "  - Loading took: %.2f seconds" %(t2 - t1)

    # If given a template for the SAM file, use it
    template = None

    if settings and "sam_template" in settings:
        template = settings["sam_template"]

    if "filter_reads" not in settings:
        filter_reads = True
    else:
        filter_reads = settings["filter_reads"]
        
    # Load the BAM file upfront
    bamfile = sam_utils.load_bam_reads(bam_filename, template=template)
        
    for gene_id, gene_info in gff_genes.iteritems():
        if gene_id not in gene_ids:
            # Skip genes that we were not asked to run on
            continue

        gene_obj = gene_info['gene_object']
        gene_hierarchy = gene_info['hierarchy']

        # Find the most inclusive transcription start and end sites for each gene
        tx_start, tx_end = gff_utils.get_inclusive_txn_bounds(gene_info['hierarchy'][gene_id])

        # Fetch reads aligning to the gene boundaries
        gene_reads = sam_utils.fetch_bam_reads_in_gene(bamfile, gene_obj.chrom,
                                                       tx_start, tx_end,
                                                       gene_obj)

        # Align the reads to the isoforms
        #reads = sam_utils.sam_reads_to_isoforms(gene_reads, gene_obj, read_len,
        #                                        overhang_len,
        #                                        paired_end=paired_end)
        reads, num_raw_reads = sam_utils.sam_parse_reads(gene_reads,
                                                         paired_end=paired_end)
                                   
        # Skip gene if none of the reads align to gene boundaries
        if filter_reads:
            if num_raw_reads < min_event_reads:
                print "Only %d reads in gene, skipping (needed >= %d reads)" \
                      %(num_raw_reads, min_event_reads)
                continue

        num_isoforms = len(gene_obj.isoforms)
        hyperparameters = ones(num_isoforms)

        ##
        ## Run the sampler
        ##
        # Create the sampler with the right parameters depending on whether
        # this is a paired-end or single-end data set.
        if paired_end:
            # Sampler parameters for paired-end mode
            sampler_params = miso.get_paired_end_sampler_params(num_isoforms,
                                                                mean_frag_len,
                                                                frag_variance,
                                                                read_len,
                                                                overhang_len=overhang_len)
            sampler = miso.MISOSampler(sampler_params, paired_end=True,
                                       log_dir=output_dir)

        else:
            # Sampler parameters for single-end mode
            sampler_params = miso.get_single_end_sampler_params(num_isoforms,
                                                                read_len,
                                                                overhang_len)
            sampler = miso.MISOSampler(sampler_params, paired_end=False,
                                       log_dir=output_dir)

        # Make directory for chromosome -- if given an event type, put
        # the gene in the event type directory
        if event_type != None:
            chrom_dir = os.path.join(output_dir, event_type, gene_obj.chrom)
        else:
            chrom_dir = os.path.join(output_dir, gene_obj.chrom)
        if not os.path.isdir(chrom_dir):
            os.makedirs(chrom_dir)
            
        output_filename = os.path.join(chrom_dir, gene_obj.label)

        sampler.run_sampler(num_iters, reads, gene_obj, hyperparameters,
                            sampler_params, output_filename, burn_in=burn_in,
                            lag=lag)
        
	    
def main():
    from optparse import OptionParser
    parser = OptionParser()
    ##
    ## Two isoform Psi
    ##
    parser.add_option("--compute-two-iso-psi", dest="two_iso_psi_files", nargs=2, default=None,
		       help="Compute Psi using MISO for a given set of two-isoform events. "
                       "Expects two arguments: the first is the set of events (in JSON/Pickle format), "
                       "the second is an output directory where estimated Psi values will "
		       "be outputted.")

    ##
    ## Multiple isoform Psi
    ##
    # parser.add_option("--compute-multi-iso-psi", dest="multi_iso_psi_files", nargs=3, default=None,
    #                   help="Compute Psi using for a given multi-isoform gene.  Expects three arguments: "
    #                   "the first is a file with the isoform lengths. The second is a file with the reads " 
    #                   "aligned to the isoform. The third is an output directory.")
    parser.add_option("--compute-gene-psi", dest="compute_gene_psi", nargs=4, default=None,
                      help="Compute Psi using for a given multi-isoform gene.  Expects four arguments: "
                      "the first is a gene ID or set of comma-separated (no spaces) gene IDs, "
                      "the second is a GFF indexed file with the gene information, the third is a sorted and "
                      "indexed BAM file with reads aligned to the gene, and the fourth is an output directory.")
    parser.add_option("--paired-end", dest="paired_end", nargs=2, default=None,
                      help="Run in paired-end mode.  Takes a mean and standard deviation "
                      "for the fragment length distribution (assumed to have discretized "
                      "normal form.)")

    ##
    ## Psi utilities
    ##
    parser.add_option("--compare-samples", dest="samples_to_compare", nargs=3, default=None,
		      help="Compute comparison statistics between the two given samples. "
                      "Expects three directories: the first is sample1's MISO output, "
                      "the second is sample2's MISO output, and the third is the directory where "
		      "results of the sample comparison will be outputted.")
    parser.add_option("--comparison-labels", dest="comparison_labels", nargs=2, default=None,
                      help="Use these labels for the sample comparison made by --compare-samples. "
                      "Takes two arguments: the label for sample 1 and the label for sample 2, "
                      "where sample 1 and sample 2 correspond to the order of samples given "
                      "to --compare-samples.")
    parser.add_option("--run-two-iso-event", dest="run_two_iso_event", nargs=3, default=None,
		      help="Run MISO on two isoform event, given an event name, an events file "
                      "(in JSON/Pickle format) and an output directory.")
    parser.add_option("--summarize-samples", dest="summarize_samples", nargs=2, default=None,
		      help="Compute summary statistics of the given set of samples. "
                      "Expects a directory with MISO output and a directory to output "
                      "summary file to.")
    parser.add_option("--summary-label", dest="summary_label", nargs=1, default=None,
                      help="Label for MISO summary file. If not given, uses basename of MISO output "
                      "directory.")
    parser.add_option("--summarize-multi-iso-samples", dest="summarize_samples", nargs=2, default=None,
		      help="Compute summary statistics of the given set of samples from multi-isoform runs. "
                      "Expects a directory with MISO output and a directory to output summary file to.")
    parser.add_option("--use-cluster", action="store_true", dest="use_cluster", default=False)
    parser.add_option("--chunk-jobs", dest="chunk_jobs", default=False, type="int",
		      help="Size (in number of events) of each job to chunk events file into. "
                      "Only applies when running on cluster.")
    parser.add_option("--settings-filename", dest="settings_filename",
                      default=os.path.join(miso_settings_path, "settings", "miso_settings.txt"),
                      help="Filename specifying MISO settings.")
    parser.add_option("--read-len", dest="read_len", type="int", default=None)
    parser.add_option("--overhang-len", dest="overhang_len", type="int", default=None)
    parser.add_option("--event-type", dest="event_type", default=None,
		      help="Event type of two-isoform events (e.g. 'SE', 'RI', 'A3SS', ...)")

    ##
    ## Gene utilities
    ##
    parser.add_option("--view-gene", dest="view_gene", nargs=1, default=None,
                      help="View the contents of a gene/event that has been indexed. "\
                      "Takes as input an indexed (.pickle) filename.")
    (options, args) = parser.parse_args()

    ##
    ## Load the settings file 
    ##
    Settings.load(os.path.expanduser(options.settings_filename))

    if options.samples_to_compare:
	sample1_dirname = os.path.abspath(options.samples_to_compare[0])
	sample2_dirname = os.path.abspath(options.samples_to_compare[1])
	output_dirname = os.path.abspath(options.samples_to_compare[2])
	if not os.path.isdir(output_dirname):
            print "Making comparisons directory: %s" %(output_dirname)
	    os.makedirs(output_dirname)
	ht.output_samples_comparison(sample1_dirname, sample2_dirname,
                                     output_dirname,
                                     sample_labels=options.comparison_labels)
	
    if options.run_two_iso_event:
	if options.read_len == None or options.overhang_len == None:
	    print "Error: must provide --read-len and --overhang-len to run."
            sys.exit(1)
            
	if options.use_cluster:
	    print "Use cluster option not supported for running on a single event."
            sys.exit(1)
            
	# convert paths to absolute path names
	event_name = options.run_two_iso_event[0]
	events_filename = os.path.abspath(options.run_two_iso_event[1]) 
	psi_outdir = os.path.abspath(os.path.expanduser(options.run_two_iso_event[2])) + '/'
        
	miso_events = as_events.MISOEvents(2, options.event_type,
                                           from_file=events_filename)
        
	run_two_iso_event(event_name, options.event_type, miso_events, psi_outdir,
			  options.read_len, options.overhang_len)

    # if options.inspect_events:
    #     print "Loading events from: %s" %(options.inspect_events)
    #     miso_events = as_events.MISOEvents(2, options.event_type, from_file=options.inspect_events)
    #     print "  - Total of %d events." %(len(miso_events.events))
	
    if options.two_iso_psi_files:
	if options.read_len == None or options.overhang_len == None:
	    print "Error: must provide --read-len and --overhang-len to run."
            sys.exit(1)

	# convert paths to absolute path names
	events_filename = os.path.abspath(options.two_iso_psi_files[0]) 
	psi_outdir = os.path.abspath(options.two_iso_psi_files[1]) + '/'
	if options.use_cluster:
	    run_two_iso_on_cluster(miso_path, events_filename, options.event_type, psi_outdir,
                                   options.read_len, options.overhang_len,
                                   chunk_jobs=options.chunk_jobs)
	else:
	    if options.chunk_jobs:
		print "Error: Chunking jobs only applies when using the --use-cluster option " \
                      "to run MISO on cluster."
                sys.exit(1)
                
	    compute_two_iso_psi(events_filename, options.event_type, psi_outdir,
				options.read_len, options.overhang_len)

    ##
    ## Multiple isoforms interface based on SAM files
    ##
    if options.compute_gene_psi != None:
        if options.read_len == None:
            print "Error: must provide --read-len."
            sys.exit(1)

        paired_end = None

        if options.paired_end != None:
            paired_end = float(options.paired_end[0]), \
                         float(options.paired_end[1])

        overhang_len = 1

        if options.overhang_len != None:
            overhang_len = options.overhang_len

        # Genes to run on from GFF
        gene_ids = options.compute_gene_psi[0].split(",")

        # GFF filename describing genes
        gff_filename = os.path.abspath(os.path.expanduser(options.compute_gene_psi[1]))

        # BAM filename with reads
        bam_filename = os.path.abspath(os.path.expanduser(options.compute_gene_psi[2]))

        # Output directory
        output_dir = os.path.abspath(os.path.expanduser(options.compute_gene_psi[3]))

        compute_gene_psi(gene_ids, gff_filename, bam_filename, output_dir,
                         options.read_len, overhang_len, paired_end=paired_end,
                         event_type=options.event_type)


    ##
    ## Summarizing samples
    ##
    if options.summarize_samples:
	samples_dir = os.path.abspath(os.path.expanduser(options.summarize_samples[0]))
        if options.summary_label != None:
            samples_label = options.summary_label
            print "Using summary label: %s" %(samples_label)
        else:
            samples_label = os.path.basename(os.path.expanduser(samples_dir))
	assert(len(samples_label) >= 1)
	summary_output_dir = os.path.abspath(os.path.join(os.path.expanduser(options.summarize_samples[1]),
							  'summary'))
	if not os.path.isdir(summary_output_dir):
	    os.makedirs(summary_output_dir)
	    
	summary_filename = os.path.join(summary_output_dir,
					'%s.miso_summary' %(samples_label))
	summarize_sampler_results(samples_dir, summary_filename)

    if options.view_gene != None:
        indexed_gene_filename = os.path.abspath(os.path.expanduser(options.view_gene))
        print "Viewing genes in %s" %(indexed_gene_filename)
        gff_genes = gff_utils.load_indexed_gff_file(indexed_gene_filename)

        if gff_genes == None:
            print "No genes."
            sys.exit(1)

        for gene_id, gene_info in gff_genes.iteritems():
            print "Gene %s" %(gene_id)
            gene_obj = gene_info['gene_object']
            print " - Gene object: ", gene_obj
            print "=="
            print "Isoforms: "
            for isoform in gene_obj.isoforms:
                print " - ", isoform
            print "=="
            print "mRNA IDs: "
            for mRNA_id in gene_info['hierarchy'][gene_id]['mRNAs']:
                print "%s" %(mRNA_id)
            print "=="    
            print "Exons: "
            for exon in gene_obj.parts:
                print " - ", exon

if __name__ == '__main__':
    main()
