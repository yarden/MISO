##
## Utilities for working with MISO output samples
##

from scipy import *
from numpy import *

import os
import glob
import misopy

from misopy.parse_csv import *
from misopy.credible_intervals import *

def maxi(l):
    m = max(l)
    for i, v in enumerate(l):
        if m == v:
            return i

def load_samples(samples_file):
    """
    Load a file with samples.  Return the samples, header from the file, the sampled MAP estimate,
    and the sampled MAP's log score.
    """
    data, h = csv2array(samples_file, skiprows=1, raw_header=True)
    sampled_map_indx = maxi(data['sampled_psi'])
    sampled_map = [float(v) for v in data['sampled_psi'][sampled_map_indx].split(',')]
    sampled_map_log_score = data['log_score'][sampled_map_indx]
#    print "  - Sampled MAP: %s" %(sampled_map)
#    print "  - Sampled MAP log score: %.4f" %(sampled_map_log_score)
    samples = []
    for vals in data['sampled_psi']:
        psi_vals = [float(v) for v in vals.split(',')]
        samples.append(psi_vals)
    samples = array(samples)
    
    # Extract counts from the file's header
    counts_info = get_counts_from_header(h[0])
    
    return (samples, h, data['log_score'], sampled_map,
            sampled_map_log_score, counts_info)


def parse_sampler_params(miso_filename):
    """
    Parse parameters that were used to produce a set of samples.  
    """
    miso_file = open(miso_filename, 'r')
    header = miso_file.readline().strip()
    miso_file.close()
    
    if header[0] == '#':
	# strip header start
	header = header[1:]
    fields = header.split('\t')
    params = {}

    for field in fields:
	key, value = field.split('=')
        params[key] = value

    return params


def get_isoforms_from_header(samples_header):
    """
    Given header of a raw MISO samples file, return the isoforms
    field.
    """
    # Get first field (removing prefix comment symbol, '#')
    isoforms = samples_header[1:].split("\t")[0]
    # Remove isoforms= key
    isoforms = isoforms.split("isoforms=")[1]
    # Remove brackets, '[', ']'
    isoforms = isoforms[1:-1]
    
    return isoforms


def get_counts_from_header(samples_header):
    """
    Given a header of a raw MISO samples file, return the
    counts= and assigned_counts= fields.
    """
    fields = samples_header[1:].split("\t")
    counts = {}
    for f in fields:
        if f.startswith("counts="):
            counts['counts'] = f.split("=")[1]
        elif f.startswith("assigned_counts="):
            counts['assigned_counts'] = f.split("=")[1]
            
    if len(counts.keys()) != 2:
        print "Warning: Could not get counts fields out of " \
              "%s header." %(samples_header)
        counts = {'counts': 'n/a',
                  'assigned_counts': 'n/a'}
        
    return counts


def summarize_sampler_results(samples_dir, summary_filename):
    """
    Given a set of samples from MISO, output a summary file.
    """
    summary_file = open(summary_filename, 'w')
    header_fields = ["event_name", "miso_posterior_mean", "ci_low", "ci_high",
                     "isoforms", "counts", "assigned_counts"]
    summary_header = "%s\n" %("\t".join(header_fields))
    summary_file.write(summary_header)
    print "Loading events from: %s" %(samples_dir)
    print "Writing summary to: %s" %(summary_filename)
    all_filenames = get_samples_dir_filenames(samples_dir)
    num_events = 0
    
    for samples_filename in all_filenames:
        basename = os.path.basename(samples_filename)
    	event_name = basename.split('.')[0]

        # Load samples and header information
	samples_results = load_samples(samples_filename)
	samples = samples_results[0]
        header = samples_results[1]
        header = header[0]

        # Get counts information from header
        counts_info = samples_results[5]
        
        num_samples, num_isoforms = shape(samples)
        output_fields = format_credible_intervals(event_name, samples)
            
        # Add isoforms information to output fields
        isoforms_field = get_isoforms_from_header(header)
        output_fields.append(isoforms_field)

        # Add counts information to output fields
        output_fields.append(counts_info['counts'])
        output_fields.append(counts_info['assigned_counts'])
        
        output_line = "%s\n" %("\t".join(output_fields))
	summary_file.write(output_line)
	num_events += 1
    print "  - Summarized a total of %d events." %(num_events)
    summary_file.close()

    
def get_samples_dir_filenames(samples_dir):
    """
    Get all the filenames associated with a samples directory.

    Assumes samples directory have the following structure:

      - samples_dir
        - chr1
        - chr2
        ...
        - chrN

    Also collect files in samples_dir for backwards compatibility.
    """
    directories = glob.glob(os.path.join(samples_dir, "*"))
    # Keep only numeric directories or directories beginning with "chr"
    directories = filter(lambda d: os.path.basename(d).startswith("chr") or \
                         os.path.basename(d).isdigit(),
                         directories)
    
    # Filenames indexed by chromosomes
    filenames = []

    for directory in directories:
        if os.path.isdir(directory):
            dir_filenames = os.listdir(directory)
            dir_filenames = [os.path.join(directory, dname) \
                             for dname in dir_filenames]
            filenames.extend(dir_filenames)

    # Filenames in top-level directory
    filenames.extend(os.listdir(samples_dir))

    # Add parent directory to all filenames
    filenames = [os.path.join(samples_dir, fname) \
                 for fname in filenames]

    # Remove directories and files beginning with "."
    filenames = filter(lambda f: not os.path.isdir(f),
                       filenames)
    filenames = filter(lambda f: not os.path.basename(f).startswith("."),
                       filenames)

    # Remove files that do not end with proper extension
    filenames = filter(lambda f: os.path.basename(f).endswith(".miso"),
                       filenames)
    return filenames


