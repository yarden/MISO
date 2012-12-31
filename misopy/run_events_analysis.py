#!/usr/bin/env python
##
## Events analysis
##

import os
import csv
import sys
import subprocess
from collections import defaultdict

import misopy.gff_utils as gff_utils
import misopy.as_events as as_events
import misopy.run_miso as run_miso
import misopy.misc_utils as misc_utils
import misopy.exon_utils as exon_utils
from misopy.parse_csv import *
from misopy.settings import Settings, load_settings
from misopy.settings import miso_path as miso_settings_path
import misopy.cluster_utils as cluster_utils

miso_path = os.path.dirname(os.path.abspath(__file__))

class GenesDispatcher:
    """
    Send MISO commands to cluster or locally
    using multi-cores.
    """
    def __init__(self, gff_dir, bam_filename,
                 output_dir, read_len, overhang_len,
                 settings_fname=None,
                 paired_end=None,
                 use_cluster=False,
                 chunk_jobs=200,
                 SGEarray=False,
                 sge_job_name="misojob"):
        self.threads = {}
        self.gff_dir = gff_dir
        self.bam_filename = bam_filename
        # Check that the BAM filename exists and that it has an index
        if not os.path.isfile(self.bam_filename):
            print "Error: BAM file %s not found." %(self.bam_filename)
            sys.exit(1)
        self.bam_index_fname = "%s.bai" %(self.bam_filename)
        if not os.path.isfile(self.bam_index_fname):
            print "WARNING: Expected BAM index file %s not found." \
                %(self.bam_index_fname)
            print "Are you sure your BAM file is indexed?"
        self.output_dir = output_dir
        self.read_len = read_len
        self.overhang_len = overhang_len
        self.settings_fname = settings_fname
        self.paired_end = paired_end
        self.use_cluster = use_cluster
        self.chunk_jobs = chunk_jobs
        self.settings = Settings.get()
        self.cluster_cmd = Settings.get_cluster_command()
        self.sge_job_name = sge_job_name
        # if chunk_jobs not given (i.e. set to False),
        # then set it to arbitrary value
        if not self.chunk_jobs:
            self.chunk_jobs = 200
        self.SGEarray = SGEarray
        self.num_processors = Settings.get_num_processors()
        self.long_thresh = 50
        self.batch_logs_dir = \
            os.path.join(output_dir, "batch-logs")
        self.batch_genes_dir = \
            os.path.join(output_dir, "batch-genes")
        self.cluster_scripts_dir = \
            os.path.join(output_dir, "cluster_scripts")
        self.scripts_output_dir = \
            os.path.join(output_dir, "scripts_output")
        misc_utils.make_dir(self.batch_logs_dir)
        misc_utils.make_dir(self.batch_genes_dir)
        misc_utils.make_dir(self.cluster_scripts_dir)
        misc_utils.make_dir(self.scripts_output_dir)
        # First compile a set of genes that should be run on
        # and output them to file along with their indexed
        # filenames
        self.gene_ids_to_gff_index = \
            gff_utils.get_gene_ids_to_gff_index(gff_dir)
        self.batch_filenames = self.output_batch_files()


    def output_batch_files(self):
        """
        Output a series of batch files containing
        gene IDs and their indexed GFF filenames.

        Return the batch filenames and their size.
        """
        batch_filenames = []
        chunk_jobs = self.chunk_jobs
        all_gene_ids = self.gene_ids_to_gff_index.keys()
        num_genes = len(all_gene_ids)
        num_chunks = max(1, int(round(num_genes / float(chunk_jobs))))
        if not self.use_cluster:
            # When not using cluster, use local multi-cores
            # using default number of processors
            num_chunks = self.num_processors
        gene_ids_batches = cluster_utils.chunk_list(all_gene_ids,
                                                    num_chunks)
        for batch_num, gene_ids_batch in enumerate(gene_ids_batches):
            batch_size = len(gene_ids_batch)
            batch_fname = os.path.join(self.batch_genes_dir,
                                       "batch-%d_genes.txt" %(batch_num))
            with open(batch_fname, "w") as batch_out:
                for gene_id in gene_ids_batch:
                    index_fname = self.gene_ids_to_gff_index[gene_id]
                    output_line = "%s\t%s\n" %(gene_id,
                                               index_fname)
                    batch_out.write(output_line)
            batch_filenames.append((batch_fname, batch_size))
        return batch_filenames


    def run(self, delay_constant=0.9):
        """
        Run batches either locally on multi-cores
        or using cluster.
        """
        batch_filenames = self.output_batch_files()
        # All MISO commands, each correspond to a batch,
        # and the number of jobs in each batch
        all_miso_cmds = []
        num_batches = len(batch_filenames)
        ##
        ## Prepare all the files necessary to run each batch
        ##
        print "Preparing to run %d batches of jobs..." %(num_batches)
        miso_run = os.path.join(miso_path, "run_miso.py")
        for batch_num, batch in enumerate(batch_filenames):
            batch_filename, batch_size = batch
            miso_cmd = \
              "python %s --compute-genes-from-file \"%s\" %s %s --read-len %d " \
                    %(miso_run,
                      batch_filename,
                      self.bam_filename,
                      self.output_dir,
                      self.read_len)
            # Add paired-end parameters and read len/overhang len
            if self.paired_end != None:
                # Run in paired-end mode
                frag_mean = float(self.paired_end[0])
                frag_sd = float(self.paired_end[1])
                miso_cmd += " --paired-end %.1f %.1f" %(frag_mean,
                                                        frag_sd)
            else:
                # Overhang len only used in single-end mode
                miso_cmd += " --overhang-len %d" %(self.overhang_len)
            # Add settings filename if given
            if self.settings_fname != None:
                miso_cmd += " --settings-filename %s" \
                    %(self.settings_fname)
            all_miso_cmds.append((miso_cmd, batch_size))
        ##
        ## Run all MISO commands for the batches
        ## either locally using multi-cores or on cluster
        ##
        # First handle special case of SGE cluster submission
        if self.use_cluster and self.SGEarray:
            print "Using SGEarray..."
            # Call SGE
            batch_argfile = os.path.join(self.cluster_scripts_dir,
                                         "run_args.txt")
            cluster_utils.run_SGEarray_cluster(all_miso_cmds,
                                               batch_argfile,
                                               self.output_dir,
                                               settings=self.settings_fname,
                                               job_name=self.sge_job_name,
                                               chunk=self.chunk_jobs)
            # End SGE case
            return
        # All cluster jobs 
        cluster_jobs = []
        for batch_num, cmd_info in enumerate(all_miso_cmds):
            miso_cmd, batch_size = cmd_info
            print "Running batch of %d genes.." %(batch_size)
            print "  - Executing: %s" %(miso_cmd)
            # Make a log file for the batch, where all the output
            # will be redirected
            time_str = time.strftime("%m-%d-%y_%H:%M:%S")
            batch_logfile = os.path.join(self.batch_logs_dir,
                                         "batch-%d-%s.log" %(batch_num,
                                                             time_str))
            cmd_to_run = "%s >> \"%s\";" %(miso_cmd, batch_logfile)
            if not self.use_cluster:
                # Run locally
                p = subprocess.Popen(cmd_to_run, shell=True)
                thread_id = "batch-%d" %(batch_num)
                print "  - Submitted thread %s" %(thread_id)
                self.threads[thread_id] = p
            else:
                # Run on cluster
                if batch_size >= self.long_thresh:
                    queue_type = "long"
                else:
                    queue_type = "short"
                # Run on cluster
                job_name = "gene_psi_batch_%d" %(batch_num)
                print "Submitting to cluster: %s" %(cmd_to_run)
                job_id = \
                    cluster_utils.run_on_cluster(cmd_to_run,
                                                 job_name,
                                                 self.output_dir,
                                                 queue_type=queue_type,
                                                 settings_fname=self.settings_fname)
                if job_id is not None:
                    cluster_jobs.append(job_id)
                time.sleep(delay_constant)
        # If ran jobs on cluster, wait for them if there are any
        # to wait on.
        cluster_utils.wait_on_jobs(cluster_jobs,
                                   self.cluster_cmd)
        # If ran jobs locally, wait on them to finish
        self.wait_on_threads()


    def wait_on_threads(self):
        if self.use_cluster:
            # If ran jobs on cluster, nothing to wait for
            return
        threads_completed = {}
        num_threads = len(self.threads)
        if num_threads == 0:
            return
        print "Waiting on %d threads..." %(num_threads)
        t_start = time.time()
        for thread_name in self.threads:
            if thread_name in threads_completed:
                continue
            curr_thread = self.threads[thread_name]
            curr_thread.wait()
            if curr_thread.returncode != 0:
                print "WARNING: Thread %s might have failed..." \
                    %(thread_name)
            threads_completed[thread_name] = True
        t_end = time.time()
        duration = ((t_end - t_start) / 60.) / 60.
        print "  - Threads completed in %.2f hours." \
            %(duration)
            

def get_ids_passing_filter(gff_index_dir,
                           bam_filename,
                           output_dir):
    """
    Apply filter to events using bedtools and return
    only the events that meet the filter.
    """
    min_reads = 20
    settings = Settings.get()
    min_event_reads = Settings.get_min_event_reads()
    
    # Check that this was indexed with a version that outputs
    # genes.gff file
    genes_gff_fname = os.path.join(gff_index_dir,
                                   "genes.gff")
    if not os.path.isfile(genes_gff_fname):
        print "WARNING: Could not find \'genes.gff\' in %s - " \
              "skipping prefilter stage. Please reindex your " \
              "GFF file with the latest version to enable " \
              "prefiltering." %(gff_index_dir)
        return None
    print "Prefiltering reads..."
    coverage_fname = exon_utils.get_bam_gff_coverage(bam_filename,
                                                     genes_gff_fname,
                                                     output_dir)
    ids_passing_filter = []
    with open(coverage_fname) as coverage_in:
        for line in coverage_in:
            # Skip comments
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            # Get the counts field and the event ID
            # if it passes the filter
            counts = int(fields[9])
            if counts < min_event_reads:
                continue
            attribs = gff_utils.parse_gff_attribs(fields[8])
            if "ID" not in attribs:
                print "WARNING: No ID= found for line:\n%s\nSkipping..." \
                    %(line)
                continue
            event_id = attribs["ID"]
            ids_passing_filter.append(event_id)
    return ids_passing_filter
            

def compute_all_genes_psi(gff_dir, bam_filename, read_len, output_dir,
                          use_cluster=False,
                          SGEarray=False,
                          chunk_jobs=800,
                          overhang_len=1,
                          paired_end=None,
                          settings_fname=None,
                          job_name="misojob",
                          prefilter=False):
    """
    Compute Psi values for genes using a GFF and a BAM filename.

    SGE functionality contributed by Michael Lovci.

    Options:
    - prefilter: if set to True, prefilter events by coverage.
      Uses bedtools to determine coverage of each event and remove
      events that do not meet the coverage criteria from the run.
    """
    print "Computing gene-level Psi for genes..." 
    print "  - GFF index: %s" %(gff_dir)
    print "  - BAM: %s" %(bam_filename)
    print "  - Read length: %d" %(read_len)
    print "  - Output directory: %s" %(output_dir)

    misc_utils.make_dir(output_dir)
    
    # Prefilter events that do not meet the coverage criteria
    # If filtering is on, only run on events that meet
    # the filter.
    if prefilter:
        print "  - Prefiltering on"
        if misc_utils.which("bedtools") is None:
            print "Error: Cannot use bedtools. Bedtools is " \
                  "required for --prefilter option"
            sys.exit(1)
        filtered_gene_ids = get_ids_passing_filter(gff_dir,
                                                   bam_filename,
                                                   output_dir)
        # Prefiltering succeeded, so process only gene ids that
        # pass the filter
        if filtered_gene_ids != None:
            num_pass = len(filtered_gene_ids)
            all_gene_ids = filtered_gene_ids
            # If none of the events meet the read coverage filter
            # something must have gone wrong, e.g. mismatch
            # in chromosome headers between BAM and GFF
            if num_pass == 0:
                print "Error: None of the events in %s appear to meet the " \
                      "read coverage filter. Check that your BAM headers " \
                      "in %s match the GFF headers of indexed events." \
                      %(gff_dir,
                        bam_filename)
                sys.exit(1)
            print "  - Total of %d events pass coverage filter." \
                %(num_pass)

    ##
    ## Submit jobs either using cluster or locally
    ## using multi-cores.
    ##
    dispatcher = GenesDispatcher(gff_dir,
                                 bam_filename,
                                 output_dir,
                                 read_len,
                                 overhang_len,
                                 settings_fname=settings_fname,
                                 paired_end=paired_end,
                                 use_cluster=use_cluster,
                                 chunk_jobs=chunk_jobs,
                                 sge_job_name=job_name,
                                 SGEarray=SGEarray)
    dispatcher.run()
            
def output_gene_ids_in_batches(gene_ids_to_gff_index,
                               output_dir,
                               chunk_jobs):
    """
    Takes mapping from gene IDs to their GFF index
    and splits them into 1 or more files (batches),
    each file containing the gene ID and its associated
    indexed GFF path.  These files are then taken as input
    by MISO and run either on a cluster or locally using
    multi-cores.
    """
    for gene_id in gene_ids_to_gff_index:
        indeed_filename = gene_ids_to_gff_index[gene_id]
        
    pass

    
        
def compute_psi(sample_filenames, output_dir, event_type,
                read_len, overhang_len,
		use_cluster=False,
                chunk_jobs=False,
                filter_events=True,
                events_info_filename=None,
                settings_filename=None):
    """
    Compute Psi values for skipped exons.  Sample filenames is a mapping from
    sample label to sample.

      - sample_filenames = [[sample_label1, sample_filename1],
                            [sample_label2, sample_filename2]]
      - output_dir: output directory
      - event_type: 'SE', 'RI', etc.
    """
    misc_utils.make_dir(output_dir)
    
    output_dir = os.path.join(output_dir, event_type)
    output_dir = os.path.abspath(output_dir)

    misc_utils.make_dir(output_dir)
	
    print "Computing Psi for events of type %s" %(event_type)
    print "  - samples used: ", sample_filenames.keys()

    for sample_label, sample_filename in sample_filenames.iteritems():
	print "Processing sample: label=%s, filename=%s" \
            %(sample_label, sample_filename)
	results_output_dir = os.path.join(output_dir, sample_label)
        misc_utils.make_dir(results_output_dir)

	# Load the set of counts and serialize them into JSON
	events = \
            as_events.load_event_counts(sample_filename,
                                        event_type,
                                        events_info_filename=events_info_filename)

	# Filter events
	if filter_events:
	    print "Filtering events..."
	    events.filter_events(settings=Settings.get())

	print "Running on a total of %d events." %(len(events.events))
	    
	events_filename = events.output_file(results_output_dir,
                                             sample_label)
	
	# Run MISO on them
	miso_cmd = "python %s --compute-two-iso-psi %s %s --event-type %s " \
                   "--read-len %d --overhang-len %d " \
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


def greeting(parser=None):
    print "MISO (Mixture of Isoforms model)"
    print "Probabilistic analysis of RNA-Seq data to detect " \
          "differential isoforms"
    print "Use --help argument to view options.\n"
    if parser is not None:
        parser.print_help()


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-genes-psi", dest="compute_genes_psi",
                      nargs=2, default=None,
                      help="Compute Psi values for a given GFF annotation "
                      "of either whole mRNA isoforms or isoforms produced by "
                      "single alternative splicing events. Expects two "
                      "arguments: an indexed GFF directory with genes to "
                      "process, and a sorted, indexed BAM file (with "
                      "headers) to run on.")
    parser.add_option("--event-type", dest="event_type", nargs=1,
		      help="Type of event (e.g. SE, RI, A3SS, ...)",
                      default=None)
    parser.add_option("--use-cluster", dest="use_cluster",
                      action="store_true", default=False,
		      help="Run events on cluster.")
    parser.add_option("--chunk-jobs", dest="chunk_jobs",
                      default=False, type="int",
		      help="Size (in number of events) of each job to chunk "
                      "events file into. Only applies when running on cluster.")
    parser.add_option("--no-filter-events", dest="no_filter_events",
                      action="store_true", default=False,
		      help="Do not filter events for computing Psi. "
		      "By default, MISO computes Psi only for events that "
                      "have a sufficient number of junction reads. "
                      "The default filter varies by event type.")
    parser.add_option("--settings-filename", dest="settings_filename",
                      default=os.path.join(miso_settings_path,
                                           "settings",
                                           "miso_settings.txt"),                    
                      help="Filename specifying MISO settings.")
    parser.add_option("--read-len", dest="read_len", default=None, type="int",
		      help="Length of sequenced reads.")
    parser.add_option("--paired-end", dest="paired_end", nargs=2, default=None, 
		      help="Run in paired-end mode. Takes mean and "
                      "standard deviation of insert length distribution.")
    parser.add_option("--overhang-len", dest="overhang_len",
                      default=None, type="int",
		      help="Length of overhang constraints "
                      "imposed on junctions.")
    parser.add_option("--output-dir", dest="output_dir", default=None,
		      help="Directory for MISO output.")
    parser.add_option("--job-name", dest="job_name", nargs=1,
                      help="Name for jobs submitted to queue for SGE jobs. " \
                      "Default is misojob", default="misojob")
    parser.add_option("--SGEarray", dest="SGEarray",
                      action="store_true", default=False,
                      help="Use MISO on cluster with Sun Grid Engine. "
                      "To be used in conjunction with --use-cluster option.")
    parser.add_option("--prefilter", dest="prefilter", default=False,
                      action="store_true",
                      help="Prefilter events based on coverage. If given as " 
                      "argument, run will begin by mapping BAM reads to event "
                      "regions (using bedtools), and omit events that do not "
                      "meet coverage criteria from the run. By default, turned "
                      "off. Note that events that do not meet the coverage criteria "
                      "will not be processed regardless, but --prefilter simply "
                      "does this filtering step at the start of the run, potentially "
                      "saving computation time so that low coverage events will not "
                      "be processed or distributed to jobs if MISO is run on a "
                      "cluster. This options requires bedtools to be installed and "
                      "available on path.")
    (options, args) = parser.parse_args()

    greeting()

    ##
    ## Load the settings file 
    ##
    if not os.path.isdir(miso_settings_path):
        print "Error: %s is not a directory containing a default MISO " \
              "settings filename. Please specify a settings filename " \
              "using --settings-filename."
        return
    
    settings_filename = \
        os.path.abspath(os.path.expanduser(options.settings_filename))
    Settings.load(settings_filename)
    
    if (not options.use_cluster) and options.chunk_jobs:
        print "Error: Chunking jobs only applies when using " \
              "the --use-cluster option to run MISO on cluster."
        sys.exit(1)
    if (not options.use_cluster) and options.SGEarray:
        print "Error: SGEarray implies that you are using an SGE cluster," \
              "please run again with --use-cluster option enabled."
        sys.exit(1)

    ##
    ## Quantitation using BAM for all genes
    ##
    if options.compute_genes_psi != None:
        # GFF filename with genes to process
        gff_filename = \
            os.path.abspath(os.path.expanduser(options.compute_genes_psi[0]))

        # BAM filename with reads
        bam_filename = \
            os.path.abspath(os.path.expanduser(options.compute_genes_psi[1]))

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
            print "WARNING: cannot use --overhang-len in paired-end mode."
            print "Using overhang = 1"

        if options.overhang_len != None:
            overhang_len = options.overhang_len
        
        compute_all_genes_psi(gff_filename, bam_filename,
                              options.read_len, output_dir,
                              overhang_len=overhang_len,
                              use_cluster=options.use_cluster,
                              SGEarray=options.SGEarray,
                              job_name=options.job_name,
                              chunk_jobs=options.chunk_jobs,
                              paired_end=options.paired_end,
                              settings_fname=settings_filename,
                              prefilter=options.prefilter)

            
		    
if __name__ == '__main__':
    main()
    
