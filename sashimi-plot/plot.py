##
## sashimi-plot
##
## Utility for visualizing RNA-Seq densities along gene models and
## for plotting MISO output
##

import os
import matplotlib

import pysam
import shelve

import GFF as gff_utils
import pe_utils

from plot_utils.samples_plotter import SamplesPlotter
from samples_utils import load_samples, parse_sampler_params
from plot_utils.plotting import *
from plot_utils.plot_gene import *
import matplotlib.pyplot as plt
from matplotlib import rc

def plot_event(event_name, pickle_dir, settings_filename,
               output_dir, png=False):
    """
    Visualize read densities across the exons and junctions
    of a given MISO alternative RNA processing event.

    Also plots MISO estimates and Psi values.
    """
    # Retrieve the full pickle filename
    genes_filename = os.path.join(pickle_dir,
                                  "genes_to_filenames.shelve")

    if not os.path.isfile(genes_filename):
        raise Exception, "Cannot find file %s. Are you sure the events " \
              "were indexed with the latest version of index_gff.py?" \
              %(genes_filename)
    
    event_to_filenames = shelve.open(genes_filename)
    if event_name not in event_to_filenames:
        raise Exception, "Event %s not found in pickled directory %s. " \
              "Are you sure this is the right directory for the event?" \
              %(event_name, pickle_dir)
    
    pickle_filename = event_to_filenames[event_name]

    # Output filename is the event name
    output_filename = os.path.join(output_dir, event_name)

    # Determine format of output file
    if png:
        format_type = "png"
    else:
        format_type = "pdf"

    output_filename += ".%s" %(format_type)
        
    print "Plotting read densities and MISO estimates along event..."
    print "  - Event: %s" %(event_name)
    print "  - Output filename: %s" %(output_filename)
        
    plot_density_from_file(pickle_filename, event_name,
                           settings_filename,
                           output_filename)


def plot_insert_len(insert_len_filename,
                    output_dir,
                    dimensions=None,
                    num_bins=25,
                    plot_sd=2,
                    png=False):
    """
    Plot insert length distribution.
    """
    if dimensions == None:
        dimensions = [7, 5]
    plt.figure(figsize=dimensions)
    if png:
        ext = "png"
    else:
        ext = "pdf"

    s = plt.subplot(1, 1, 1)
    plot_name = os.path.basename(insert_len_filename)
    output_filename = os.path.join(output_dir,
                                   "%s.%s" %(plot_name,
                                             ext))
    print "Plotting insert length distribution..."
    print "  - Distribution file: %s" %(insert_len_filename)
    print "  - Output plot: %s" %(output_filename)
    insert_dist, params = pe_utils.load_insert_len(insert_len_filename)

    mean, sdev, dispersion, num_pairs \
          = pe_utils.compute_insert_len_stats(insert_dist)
    print "min insert: %.1f" %(min(insert_dist))
    print "max insert: %.1f" %(max(insert_dist))
    plt.title("%s (%d read-pairs)" \
              %(plot_name,
                num_pairs),
              fontsize=10)
    plt.hist(insert_dist, bins=num_bins, color='k',
             edgecolor="#ffffff", align='mid')
    axes_square(s)
    ymin, ymax = s.get_ylim()
    plt.text(0.05, 0.95, "$\mu$: %.1f\n$\sigma$: %.1f\n$d$: %.1f" \
             %(round(mean, 2),
               round(sdev, 2),
               round(dispersion, 2)),
             horizontalalignment='left',
             verticalalignment='top',
             bbox=dict(edgecolor='k', facecolor="#ffffff",
                       alpha=0.5),
             fontsize=10,
             transform=s.transAxes)
    plt.xlabel("Insert length (nt)")
    plt.ylabel("No. read pairs")
    plt.savefig(output_filename)
        

def plot_posterior(miso_filename, output_dir,
                   with_intervals=None,
                   dimensions=None,
                   plot_mean=False,
                   png=False):
    """
    Plot posterior distribution.
    """
#    samples, log_scores, params, gene = load_samples(miso_filename)
    samples, h, log_scores, sampled_map,\
             sampled_map_log_score, counts_info = load_samples(miso_filename)
    params = parse_sampler_params(miso_filename)
    
    sp = SamplesPlotter(samples, params)
    
    if with_intervals != None:
        with_intervals = float(with_intervals)/100.
        print "Plotting with %d-percent confidence intervals" %(int(with_intervals * 100))
    else:
        with_intervals = False

    if plot_mean:
        print "Plotting mean of posterior."

    print "Plotting posterior distribution..."
    print "  - MISO event file: %s" %(miso_filename)
    print "  - Output dir: %s" %(output_dir)
    
    sp.plot(plot_intervals=with_intervals, fig_dims=dimensions,
            plot_mean=plot_mean)

    # Determine output format type
    if not png:
        matplotlib.use('PDF')
        plt.rcParams['ps.useafm'] = True
        plt.rcParams['pdf.fonttype'] = 42
        file_ext = ".pdf"
    else:
        file_ext = ".png"

    output_filename = os.path.join(output_dir,
                                   os.path.basename(miso_filename).replace(".miso",
                                                                           file_ext))
    print "Outputting plot to: %s" %(output_filename)
    plt.savefig(output_filename)
    

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--plot-posterior", dest="plot_posterior", nargs=1, default=None,
                      help="Plot the posterior distribution. Takes as input a raw MISO output "
                      "file (.miso)")
    parser.add_option("--plot-insert-len", dest="plot_insert_len", nargs=1, default=None,
                      help="Plot the insert length distribution from a given insert length (*.insert_len) "
                      "filename.")
    parser.add_option("--plot-event", dest="plot_event", nargs=3, default=None,
                      help="Plot read densities and MISO inferences for a given alternative event. "
                      "Takes the arguments: (1) event name (i.e. the ID= of the event based on MISO gff3 "
                      "annotation file, (2) directory where MISO output is for that event type (e.g. if event is a "
                      "skipped exon, provide the directory where the output for all SE events are), "
                      "(3) path to plotting settings file.")
    parser.add_option("--with-intervals", dest="with_intervals", nargs=1, default=None, 
                      help="Include confidence intervals in plot. To be used with --plot-posterior. "
                      "Takes an argument between 1 and 99 corresponding to the confidence "
                      "interval to be used, e.g.: 95")
    parser.add_option("--plot-mean", dest="plot_mean", action="store_true", default=False,
                      help="Plot the mean of the posterior distribution. To be used with --plot-posterior.")
    parser.add_option("--dimensions", dest="dimensions", nargs=2, default=None,
                      help="Dimensions of the outputted figure: takes width by height (in inches).")
    parser.add_option("--png", dest="png", default=False, action="store_true",
                      help="Output plot in PNG format (the default is PDF).")    
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.output_dir == None:
        print "Error: need --output-dir"
        return

    output_dir = os.path.abspath(os.path.expanduser(options.output_dir))
    dimensions = None

    if options.dimensions != None:
        dimensions = (int(options.dimensions[0]),
                      int(options.dimensions[1]))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.plot_insert_len != None:
        insert_len_filename = os.path.abspath(os.path.expanduser(options.plot_insert_len))
        plot_insert_len(insert_len_filename, output_dir, png=options.png)

    if options.plot_posterior != None:
        miso_filename = os.path.abspath(os.path.expanduser(options.plot_posterior))

        with_intervals = None
        plot_mean = options.plot_mean
        
        if options.with_intervals != None:
            with_intervals = options.with_intervals

        plot_posterior(miso_filename, output_dir,
                       with_intervals=with_intervals,
                       dimensions=dimensions,
                       plot_mean=plot_mean,
                       png=options.png)

    if options.plot_event != None:
        # Plot a MISO event along with its RNA-Seq read densities
        event_name = options.plot_event[0]
        pickle_dir = os.path.abspath(os.path.expanduser(options.plot_event[1]))
        settings_filename = os.path.abspath(os.path.expanduser(options.plot_event[2]))
        plot_event(event_name, pickle_dir, settings_filename, output_dir,
                   png=options.png)
        

if __name__ == '__main__':
    main()
