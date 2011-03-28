##
## Plotting utils 
##
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import SubplotZero
import scipy.stats as stats
from scipy import *
from scipy import linalg
import sys

def show_spines(ax,spines):
    for loc, spine in ax.spines.iteritems():
        if loc not in spines:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

def fit_line(x, y, plot_line=False):
    """
    Plot best fit least squares line, return R^2
    """
    A = vstack([x, ones(len(x))]).T
    m, c = linalg.lstsq(A, y)[0]
    if plot_line:
        plt.plot(x, m*x + c, 'r', lw=1.2)
    return square(stats.pearsonr(x, y)), m, c

def remove_extra_ticks(ax):
    for i, line in enumerate(ax.get_xticklines() + ax.get_yticklines()):
        if i%2 == 1:   # odd indices
            line.set_visible(False)

def axes_square(plot_handle):
    plot_handle.axes.set_aspect(1/plot_handle.axes.get_data_ratio()) 

def setup_two_axes(fig, labelpad=1, invisible=["bottom", "top", "right"]):
    plt.rcParams['xtick.major.pad'] = 0.1
    plt.rcParams['xtick.minor.pad'] = 0.1
    plt.rcParams['ytick.major.pad'] = 2
    plt.rcParams['ytick.minor.pad'] = 2
    ax = SubplotZero(fig, 1, 1, 1)
    ax.yaxis.labelpad = labelpad
    fig.add_subplot(ax)
    # make xzero axis (horizontal axis line through y=0) visible.
    ax.axis["xzero"].set_visible(True)
    ax.xaxis.labelpad = labelpad    
    # make other axis (bottom, top, right) invisible.
    for n in invisible:
        ax.axis[n].set_visible(False)
    return ax

def setup_two_axes_subplot(fig, m, n, curr_plot_num, invisible=["bottom", "top", "right"]):
    ax = SubplotZero(fig, m, n, curr_plot_num)
    fig.add_subplot(ax)
    ax.axis["xzero"].set_visible(True)
    for n in invisible:
        ax.axis[n].set_visible(False)
    return ax

def restyle_ticks(c, min_val=0, max_val=1):
    plt.xlim(min_val - c, max_val + c)
    plt.ylim(min_val - c, max_val + c)

def label_stacked_bars(rects1, rects2, labels, h=1.02):
    label_ind = 0
    for rect, rect_other in zip(rects1, rects2):
        height = rect.get_height() + rect_other.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., h*height, labels[label_ind], ha='center', va='bottom')
        label_ind += 1

import matplotlib.transforms as mtransforms

def expand_subplot(ax, num2):
    update_params_orig = ax.update_params

    ax._num2 = num2 - 1
    def _f(self=ax):
        num_orig = self._num

        self._num = self._num2
        update_params_orig()
        right, top = self.figbox.extents[2:]

        self._num = num_orig
        update_params_orig()
        left, bottom = self.figbox.extents[:2]

        self.figbox = mtransforms.Bbox.from_extents(left, bottom,
                                                    right, top)

    ax.update_params = _f
    ax.update_params()
    ax.set_position(ax.figbox)

colors = {'steelblue': '#63B8FF',
          'lightblue1': '#3399FF',
          'signblue': '#003F87', # darkblue
          'grey1': '#999999',
          'darkred': '#990000'}
