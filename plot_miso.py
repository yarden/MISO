##
## Utility for visualizing the output of MISO.
##

import os
import matplotlib

import pysam
import shelve
import GFF as gff_utils
from plot_utils.samples_plotter import SamplesPlotter, load_samples
from plot_utils.plotting import *
import matplotlib.pyplot as plt
from matplotlib import rc
from plot_gene import *


def plot_event(event_name, pickle_dir, settings_filename,
               output_dir, png=False):
    """
    Visualize read densities across the exons and junctions
    of a given MISO alternative RNA processing event.

    Also plots MISO estimates and Psi values.
    """
    # Retrieve the full pickle filename
    event_to_filenames = shelve.load(os.path.join(pickle_dir,
                                                  "genes_to_filenames.shelve"))
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
        
    plot_density_from_file(pickle_filename, event, settings_filename,
                           output_filename)


def plot_posterior(miso_filename, output_dir,
                   with_intervals=None,
                   dimensions=None,
                   plot_mean=False,
                   png=False):
    """
    Plot posterior distribution.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    samples, log_scores, params, gene = load_samples(miso_filename)
    sp = SamplesPlotter(samples, gene, params)
    
    if with_intervals != None:
        with_intervals = float(with_intervals)/10.
        print "Plotting with %.2f confidence intervals" %(with_intervals)
    else:
        with_intervals = False

    if plot_mean:
        print "Plotting mean of posterior."

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
    parser.add_option("--plot-event", dest="plot_event", nargs=2, default=None,
                      help="Plot read densities and MISO inferences for a given alternative event. "
                      "Takes the arguments: (1) event name (i.e. the ID= of the event based on MISO gff3 "
                      "annotation file, (2) directory where the pickled annotation files are, "
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

    if options.plot_posterior != None:
        miso_filename = os.path.abspath(os.path.expanduser(options.plot_posterior))

        with_intervals = None
        dimensions = None
        plot_mean = options.plot_mean
        
        if options.with_intervals != None:
            with_intervals = options.with_intervals[0]

        if options.dimensions != None:
            dimensions = (int(options.dimensions[0]),
                          int(options.dimensions[1]))

        plot_posterior(miso_filename, output_dir,
                       with_intervals=with_intervals,
                       dimensions=dimensions,
                       plot_mean=plot_mean,
                       png=options.png)

        pickle_filename = os.path.abspath(os.path.expanduser(options.plot_density[0]))
        bam_filename = os.path.abspath(os.path.expanduser(options.plot_density[1]))
        plot_density(pickle_filename, bam_filename)

    if options.plot_event != None:
        # Plot a MISO event along with its RNA-Seq read densities
        event_name = options.plot_event[0]
        pickle_dir = os.path.abspath(os.path.expanduser(options.plot_event[1]))
        settings_filename = os.path.abspath(os.path.expanduser(options.plot_event[2]))
        plot_event(event_name, pickle_dir, settings_filename, output_dir,
                   png=options.png)
        

if __name__ == '__main__':
    main()
