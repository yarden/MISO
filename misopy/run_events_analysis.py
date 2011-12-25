#!/usr/bin/env python
##
## Events analysis
##

import os
import csv
import sys
from collections import defaultdict

import misopy.gff_utils as gff_utils
import misopy.as_events as as_events
import misopy.run_miso as run_miso
from misopy.parse_csv import *
from misopy.settings import Settings, load_settings
import misopy.cluster_utils as cluster_utils

miso_path = os.path.dirname(os.path.abspath(__file__))

def compute_all_genes_psi(gff_dir, bam_filename, read_len, output_dir,
                          use_cluster=False, SGEarray=False, chunk_jobs=200,
                          overhang_len=1, paired_end=None,
                          settings=None, job_name="misojob"):
    """
    Compute Psi values for genes using a GFF and a BAM filename.

    SGE functionality contributed by Michael Lovci.
    """
    gene_ids_to_gff_index = gff_utils.get_gene_ids_to_gff_index(gff_dir)

    num_genes = len(gene_ids_to_gff_index.keys())

    miso_run = os.path.join(miso_path, "run_miso.py")

    print "Computing gene-level Psi for %d genes..." \
          %(num_genes)
    print "  - GFF index: %s" %(gff_dir)
    print "  - BAM: %s" %(bam_filename)
    print "  - Read length: %d" %(read_len)
    print "  - Output directory: %s" %(output_dir)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    # All commands to run
    all_miso_cmds = []

    for gene_id, gff_index_filename in gene_ids_to_gff_index.iteritems():
        miso_cmd = "python %s --compute-gene-psi \"%s\" %s %s %s --read-len %d " \
                   %(miso_run, gene_id, gff_index_filename, bam_filename, output_dir,
                     read_len)
        
        if paired_end != None:
            # Run in paired-end mode
            frag_mean = float(paired_end[0])
            frag_sd = float(paired_end[1])
            miso_cmd += " --paired-end %.1f %.1f" %(frag_mean, frag_sd)
        else:
            miso_cmd += " --overhang-len %d" %(overhang_len)

        # Add settings filename if given
        if settings != None:
            miso_cmd += " --settings-filename %s" %(settings)

        if use_cluster:
            # If asked to use cluster, accumulate the MISO commands
            # but do not run them
            all_miso_cmds.append(miso_cmd)
        else:
            print "  - Executing: %s" %(miso_cmd)
            os.system(miso_cmd)

    miso_settings = Settings.load(settings)

    if use_cluster:
        if SGEarray:
            if not chunk_jobs:
                chunk_jobs = 2500
                print "  - Using default chunk jobs = %d" %(chunk_jobs)
            cluster_output_dir = os.path.join(output_dir, "cluster_scripts")
            if not os.path.isdir(cluster_output_dir):
                os.makedirs(cluster_output_dir)
            batch_argfile = os.path.join(cluster_output_dir, "run_args.txt")
            # Call SGE if asked
            cluster_utils.run_SGEarray_cluster(all_miso_cmds, batch_argfile,
                                               output_dir,
                                               settings=settings,
                                               job_name=job_name,
                                               chunk=chunk_jobs)
        else:
            # Threshold for putting jobs in the long queue
            long_thresh = 50
            
            # Delay between jobs
            delay_constant = 0.9
        
            # Invoke the commands using the cluster
            print "Sending %d genes to be run on cluster in chunks of %d..." \
                %(num_genes, chunk_jobs)

            if not chunk_jobs:
                print "  - Using default chunk jobs = %d" %(200)
                chunk_jobs = 200

            chunk_jobs = max(1, int(round(num_genes / float(chunk_jobs))))

            # Split the gene records into batches
            cmd_batches = cluster_utils.chunk_list(all_miso_cmds, chunk_jobs)

            time_str = time.strftime("%m-%d-%y_%H:%M:%S")

            for batch_num, batch in enumerate(cmd_batches):
                batch_size = len(batch)
                print "Running batch %d (batch size = %d)" %(batch_num,
                                                             batch_size)

                if batch_size >= long_thresh:
                    queue_type = "long"
                else:
                    queue_type = "short"

                # Pool all the MISO commands belonging to this batch
                batch_logs_dir = os.path.join(output_dir, "batch-logs")
                if not os.path.isdir(batch_logs_dir):
                    os.makedirs(batch_logs_dir)
                batch_logfile = os.path.join(batch_logs_dir,
                                             "batch-%d-%s.log" %(batch_num,
                                                                 time_str))
                redirected_output = " >> %s;\n" %(batch_logfile)
                cmd_to_run = redirected_output.join(batch)

                # Run on cluster
                job_name = "gene_psi_batch_%d" %(batch_num)
                cluster_utils.run_on_cluster(cmd_to_run, job_name, output_dir,
                                             queue_type=queue_type,
                                             settings=settings)
                time.sleep(delay_constant)
        
            
        
def compute_psi(sample_filenames, output_dir, event_type, read_len, overhang_len,
		use_cluster=False, chunk_jobs=False, filter_events=True,
                events_info_filename=None, settings_filename=None):
    """
    Compute Psi values for skipped exons.  Sample filenames is a mapping from
    sample label to sample.

      - sample_filenames = [[sample_label1, sample_filename1],
                            [sample_label2, sample_filename2]]
      - output_dir: output directory
      - event_type: 'SE', 'RI', etc.
    """
    if not os.path.isdir(output_dir):
	os.makedirs(output_dir)
	
    output_dir = os.path.join(output_dir, event_type)
    output_dir = os.path.abspath(output_dir)
    if not os.path.isdir(output_dir):
	os.makedirs(output_dir)
	
    print "Computing Psi for events of type %s" %(event_type)
    print "  - samples used: ", sample_filenames.keys()
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    for sample_label, sample_filename in sample_filenames.iteritems():
	print "Processing sample: label=%s, filename=%s" %(sample_label, sample_filename)
	results_output_dir = os.path.join(output_dir, sample_label)
        if not os.path.isdir(results_output_dir):
            os.makedirs(results_output_dir)

	# Load the set of counts and serialize them into JSON
	events = as_events.load_event_counts(sample_filename, event_type,
                                             events_info_filename=events_info_filename)

	# Filter events
	if filter_events:
	    print "Filtering events..."
	    events.filter_events(settings=Settings.get())

	print "Running on a total of %d events." %(len(events.events))
	    
	events_filename = events.output_file(results_output_dir, sample_label)
	
	# Run MISO on them
	miso_cmd = 'python %s --compute-two-iso-psi %s %s --event-type %s --read-len %d --overhang-len %d ' \
                   %(os.path.join(miso_path, 'run_miso.py'),
                     events_filename,
                     results_output_dir,
                     event_type,
                     read_len,
                     overhang_len)
	if use_cluster:
	    if chunk_jobs:
		miso_cmd += ' --use-cluster --chunk-jobs %d' %(chunk_jobs)
	    else:
		miso_cmd += ' --use-cluster'
        print "Executing: %s" %(miso_cmd)
	if use_cluster:
	    print " - Using cluster"
	os.system(miso_cmd)


def pool_counts(counts_filenames, event_type, output_dir,
                counts_delimiter=';', column_delimiter='\t',
		NA_delimiter='n/a'):
    """
    Pool together the counts in counts_filenames, save output to output_dir.
    """
    print "Pooling together %d count filenames into: %s" %(len(counts_filenames),
                                                           output_dir)
    print "  - Events type: %s" %(event_type)
    
    if not os.path.isdir(output_dir):
        print "Creating directory: %s" %(output_dir)
        os.makedirs(output_dir)

    # Compute labels for each counts file
    labels_to_filenames = {}

    all_sample_labels = []
    
    for counts_filename in counts_filenames:
        counts_filename = os.path.abspath(os.path.expanduser(counts_filename))
        if not counts_filename.endswith(counts_filename):
            print "Skipping file %s, it does not end with .counts" %(counts_filename)
            continue

        counts_label = os.path.basename(counts_filename).split('.counts')
        
        if len(counts_label) != 2:
            print "Skipping file %s, does not look like counts file" %(counts_filename)
            continue

        counts_label = counts_label[0]
        all_sample_labels.append(counts_label)
        labels_to_filenames[counts_label] = counts_filename

    # Set of ordered labels which will determine the order
    # in which counts are outputted
    fieldnames = ['event_name']
    fieldnames.extend(all_sample_labels)
    fieldnames.append('sum')
    
    pooled_counts_data = defaultdict(dict)
    
    for counts_label, counts_filename in labels_to_filenames.iteritems():
        counts_file = open(counts_filename)
        counts_data = csv.reader(counts_file, delimiter=column_delimiter)

        for event_name, event_counts in counts_data:
            counts = [int(c) for c in event_counts.split(counts_delimiter)]
            if len(counts) <= 1:
                raise Exception, "Malformed counts file: %s in %s" %(event_counts,
                                                                     counts_filename)
            pooled_counts_data[event_name][counts_label] = counts
        counts_file.close()

    # Serialize counts data
    pooled_counts_filename = os.path.join(output_dir,
                                          "%s.pooled_counts" %("_".join(all_sample_labels)))
    pooled_counts_file = open(pooled_counts_filename, 'w')
    
    # Write header
    header = column_delimiter.join(fieldnames)
    pooled_counts_file.write(header + "\n")

    pooled_counts_out = csv.writer(pooled_counts_file, delimiter=column_delimiter,
                                   quoting=csv.QUOTE_NONE)

    # Compile all counts for the event as strings; if a count
    # is unavailable, use NA_delimiter to denote it
    for event_name, event_counts in pooled_counts_data.iteritems():
        event_values = [event_name]
        all_sample_counts = []
        for sample_label in all_sample_labels:
            if sample_label not in event_counts:
                event_values.append(NA_delimiter)
            else:
                # Compile all of the event's counts from the sample
                curr_sample_counts = [c for c in event_counts[sample_label]]
                all_sample_counts.append(curr_sample_counts)

                # Convert them to strings for outputting to file
                curr_sample_counts_str = [str(c) for c in curr_sample_counts]
                event_values.append(counts_delimiter.join(curr_sample_counts_str))

        # Sum all the event's counts from the samples in which it is non-NA
        sum_of_counts = [str(c) for c in sum(all_sample_counts, 0)]

        # Format the summed counts
        event_values.append(counts_delimiter.join(sum_of_counts))
        pooled_counts_out.writerow(event_values)
    return pooled_counts_data


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-events-psi", dest="compute_events_psi", nargs=2, default=None,
		      help="Compute Psi values for all events. Expects two arguments: a "
		      "set of labels and a set of filenames with associated "
		      "read counts.")
    parser.add_option("--compute-genes-psi", dest="compute_genes_psi", nargs=2, default=None,
                      help="Compute Psi values for a given GFF annotation of either whole mRNA isoforms "
                      "or isoforms produced by single alternative splicing events. "
                      "Expects two arguments: an indexed GFF directory with genes to process, and a sorted, indexed "
                      "BAM file (with headers) to run on.")
    parser.add_option("--event-type", dest="event_type", nargs=1,
		      help="Type of event (e.g. SE, RI, A3SS, ...)",
                      default=None)
    parser.add_option("--events-file", dest="events_info_filename", nargs=1,
		      help="Filename with all the events and their coordinates. This is used "
		      "to compute the length of the alternative regions in each event.",
		      default=None)
    parser.add_option("--use-cluster", dest="use_cluster", action="store_true", default=False,
		      help="Run events on cluster.")
    parser.add_option("--chunk-jobs", dest="chunk_jobs", default=False, type="int",
		      help="Size (in number of events) of each job to chunk events file into. "
		      "Only applies when running on cluster.")
    parser.add_option("--no-filter-events", dest="no_filter_events", action="store_true",
		      default=False,
		      help="Do not filter events for computing Psi. "
		      "By default, MISO computes Psi only for events that have a "
		      "sufficient number of junction reads. The default filter "
		      "varies by event type.")
    parser.add_option("--settings-filename", dest="settings_filename",
                      default=os.path.join(miso_path, "settings", "miso_settings.txt"),                    
                      help="Filename specifying MISO settings.")
    parser.add_option("--pool-counts", dest="pool_counts", default=None,
		      help="Given a series of comma separated filenames (no spaces), "
		      "pool them together. Must specify an output directory "
		      "with --output-dir.")
    parser.add_option("--read-len", dest="read_len", default=None, type="int",
		      help="Length of sequenced reads.")
    parser.add_option("--paired-end", dest="paired_end", nargs=2, default=None, 
		      help="Run in paired-end mode. Takes mean and standard deviation "
                      "of insert length distribution.")
    parser.add_option("--overhang-len", dest="overhang_len", default=None, type="int",
		      help="Length of overhang constraints imposed on junctions.")
    parser.add_option("--output-dir", dest="output_dir", default=None,
		      help="Directory for MISO output.")
    parser.add_option("--job-name", dest="job_name", nargs=1, help="name for jobs submitted to queue. default is misojob", default="misojob")
    parser.add_option("--SGEarray", dest="SGEarray", action="store_true", default=False)
    (options, args) = parser.parse_args()

    ##
    ## Load the settings file 
    ##
    settings_filename = os.path.abspath(os.path.expanduser(options.settings_filename))
    Settings.load(settings_filename)
    
    print "Loading settings file from: %s" %(settings_filename)


    if (not options.use_cluster) and options.chunk_jobs:
        print "Error: Chunking jobs only applies when using " \
              "the --use-cluster option to run MISO on cluster."
        sys.exit(1)
    if (not options.use_cluster) and options.SGEarray:
        print "Error: SGEarray implies that you are using an SGE cluster," \
              "please run again with --use-cluster option enabled."
        sys.exit(1)
    
    if options.pool_counts:
        if options.output_dir == None:
            print "Error: Need output directory to pool counts together. Use --output-dir."
            sys.exit(1)

        if options.event_type == None:
            print "Error: Need event type to perform pooling."
            sys.exit(1)

        
        counts_filenames = [os.path.normpath(os.path.expanduser(count_filename)) \
                            for count_filename in options.pool_counts.split(',')]

        if len(counts_filenames) < 2:
            print "Error: Need 2 filenames or more to perform pooling."
            sys.exit(1)
        
        print "Pooling together: "
        for filename in counts_filenames:
            print "  - %s" %(filename)
            
        print "Outputting pooled counts to: %s" %(options.output_dir)
        pool_counts(counts_filenames, options.event_type,
                    os.path.expanduser(options.output_dir))

    ##
    ## Event types that require additional event files
    ##
    event_types_requiring_files = ['AFE', 'ALE']

    if options.compute_events_psi:
        print "Computing Psi for events..."
	# Error check arguments
	if not options.event_type:
	    print "Error: Need event type (e.g. SE, TandemUTR, ...) to run MISO."
            sys.exit(1)

	if (options.event_type in event_types_requiring_files) and \
	   options.events_info_filename == None:
	    print "Error: Need for event information filename for events of type %s" \
		  %(options.event_type)
            sys.exit(1)
	
	if not (options.read_len != None and options.overhang_len != None):
	    print "Error: Need read length and overhang length to run MISO."
            sys.exit(1)

	if options.output_dir == None:
	    print "Error: Need output directory to run MISO."
            sys.exit(1)


	labels = options.compute_events_psi[0].split(',')
	filenames = options.compute_events_psi[1].split(',')
	assert(len(labels) == len(filenames))
	
	sample_filenames = {}
	for label, filename in zip(labels, filenames):
	    sample_filenames[label] = os.path.normpath(os.path.expanduser(filename))

	filter_events = not options.no_filter_events

        print "Filter events? "
        if filter_events:
            print "  - yes"
        else:
            print "  - no"

        events_info_filename = None
        if options.events_info_filename != None:
            events_info_filename = os.path.expanduser(options.events_info_filename)

	compute_psi(sample_filenames, os.path.expanduser(options.output_dir), options.event_type,
                    options.read_len, options.overhang_len, use_cluster=options.use_cluster,
                    filter_events=filter_events, chunk_jobs=options.chunk_jobs,
                    events_info_filename=events_info_filename,
                    settings_filename=settings_filename)

    ##
    ## Gene-level isoform quantitation using BAM
    ##
    if options.compute_genes_psi != None:
        # GFF filename with genes to process
        gff_filename = os.path.abspath(os.path.expanduser(options.compute_genes_psi[0]))

        # BAM filename with reads
        bam_filename = os.path.abspath(os.path.expanduser(options.compute_genes_psi[1]))

        if options.output_dir == None:
            print "Error: need --output-dir to compute gene-level Psi."
            sys.exit(1)

        # Output directory to use
        output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

        if options.read_len == None:
            print "Error: need --read-len when computing gene-level Psi."
            sys.exit(1)

        overhang_len = 1

        if options.paired_end != None and options.overhang_len != None:
            print "Error: cannot use --overhang-len in paired-end mode."
            sys.exit(1)

        if options.overhang_len != None:
            overhang_len = options.overhang_len
        
        compute_all_genes_psi(gff_filename, bam_filename, options.read_len, output_dir,
                              overhang_len=overhang_len,
                              use_cluster=options.use_cluster,
                              SGEarray=options.SGEarray,
                              job_name=options.job_name,
                              chunk_jobs=options.chunk_jobs,
                              paired_end=options.paired_end,
                              settings=settings_filename)
            
		    
if __name__ == '__main__':
    main()
    
