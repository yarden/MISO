##
## Class for representing figures
##
import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

def Sashimi:
    """
    Representation of a figure.
    """
    def __init__(self, label, output_dir, dimensions=None, png=False,
                 output_filename=None, settings=None):
        """
        Initialize image settings.
        """
        # Plot label, will be used in creating the plot
        # output filename
        self.label = label

        # Set output directory
        self.set_output_dir(output_dir)

        # Plot settings
        self.settings = settings

        if output_filename != None:
            # If explicit output filename is given to us, use it
            self.output_filename = output_filename
        else:
            # Otherwise, use the label and the output directory
            self.set_output_filename()
        
        if dimensions != None:
            self.dimensions = dimensions
        else:
            self.dimensions = [7, 5]

        if png:
            self.output_ext = ".png"
        else:
            self.output_ext = ".pdf"

    def set_output_dir(self, output_dir):
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))

    def set_output_filename(self):
        plot_basename = "%s%s" %(label, self.output_ext)
        self.output_filename = os.path.join(self.output_dir, plot_basename)

    def setup_figure(self):
        plt.figure(figsize=self.dimensions)

    def save_plot(self):
        """
        Save plot to the output directory. Determine
        the file type.
        """
        if self.output_filename == None:
            raise Exception, "sashimi_plot does not know where to save the plot."
        plt.savefig(self.output_filename)

        
        
