##
## Visualize samples produced by MISO.
##

from scipy import *
from numpy import *

import matplotlib
from plot_utils.plotting import *
import matplotlib.pyplot as plt
from matplotlib import rc
from parse_csv import csv2array
import time

import Gene
import hypothesis_test as ht

class SamplesPlotter:
    """
    Visualize a set of samples from a run of MISO.
    """
    def __init__(self, samples, gene, params, log_scores=None,
                 percent_acceptance=None):
	"""
	Given a sampler instance, store its properties.
	"""
	# sampler parameters
	self.samples = samples	
	self.params = params
	self.gene = gene
	self.log_scores = log_scores
	self.percent_acceptance = percent_acceptance
	
	assert(len(samples) > 1)	
	if gene:
	    assert(len(gene.isoforms) > 1)

    def plot(self, fig=None, output_dir=None, num_rows=1, num_cols=1, subplot_start=1, title=None,
	     plot_intervals=None, value_to_label=None, label=None, bins=10, bbox_coords=None, vanilla=False,
             plot_mean=False, fig_dims=(5, 5)):
	"""
	Plot a set of samples.

	 - credible_intervals: if set to true, plot Bayesian confidence intervals 
	"""
	plot_handle = None
	if len(self.gene.isoforms) == 2:
	    plot_handle = self.plot_two_iso_samples(fig=fig, plots_dir=output_dir, num_cols=num_cols,
						    num_rows=num_rows, subplot_start=subplot_start,
						    plot_intervals=plot_intervals,
						    value_to_label=value_to_label,
						    label=label, bbox_coords=bbox_coords, title=title, vanilla=vanilla,
                                                    plot_mean=plot_mean, fig_dims=fig_dims)
	else:
	    num_isoforms = self.samples.shape[1] 
	    num_rows = 1
	    num_cols = num_isoforms
	    for c in range(num_cols):
		plot_handle = self.plot_two_iso_samples(fig, isoform_index=c, subplot_start=c + 1, num_cols=num_cols,
                                                        plot_intervals=plot_intervals,
							title=title, bins=bins, vanilla=vanilla,
                                                        plot_mean=plot_mean, fig_dims=fig_dims)
		plt.ylabel('Frequency (Isoform %d)' %(c + 1))
	    plt.subplots_adjust(wspace=0.5)
	return plot_handle
	
    def plot_two_iso_samples(self, fig=None, isoform_index=0, num_rows=1, num_cols=1, subplot_start=1,
			     plots_dir=None, map_estimate=None, simulation_num=1,
			     plot_intervals=False, value_to_label=None, label=None, plot_filename=None,
                             bins=None, bbox_coords=None, with_legend=True, title=None, vanilla=False,
                             plot_mean=False, normed=False, fig_dims=(5, 5)):
	"""
	Plot a set of samples for Psi of a two isoform gene.
	"""
	if not fig:
	    sampled_psi_fig = plt.figure(figsize=fig_dims, dpi=300)
	else:
	    sampled_psi_fig = fig
	ax = sampled_psi_fig.add_subplot(num_rows, num_cols, subplot_start)
	num_iters = int(self.params['iters'])
	burn_in = int(self.params['burn_in'])
	lag = int(self.params['lag'])
	percent_acceptance = float(self.params['percent_accept'])
	proposal_type = self.params['proposal_type']
	plt.rcParams['font.size'] = 10
	show_spines(ax, ['left', 'bottom'])
	bins = bins
	assert((value_to_label == None and label == None) or \
	       (value_to_label != None and label != None))
	# retrieve samples
	samples_to_plot = self.samples[:, isoform_index]
	# picasso blue #0276FD
        
	if not vanilla:
	    if bins != None:
		plt.hist(samples_to_plot, align='center', lw=0.5, facecolor='#0276FD',
                         edgecolor='#ffffff')
	    else:
		plt.hist(samples_to_plot, align='center', lw=0.5, facecolor='#0276FD',
                         edgecolor='#ffffff')
	else:
	    plt.hist(samples_to_plot, align='center', facecolor='#0276FD', edgecolor='#0276FD')
	plt.xlabel(r'${\hat{\Psi}}_{\mathregular{MISO}}$')
	plt.ylabel('Frequency')
	plt.xlim([0, 1])

        # Normalize samples
        if normed:
            yticks = list(plt.gca().get_yticks())
            print "yticks: ", yticks
            ytick_labels = ["%.2f" %(float(ytick) / float(normed)) for ytick in yticks]
            ax.set_yticklabels(ytick_labels)
            
	curr_axes = plt.gca()
	# Plot MAP estimate for same data
	if map_estimate:
	    l = plt.axvline(x=map_estimate, color='b', linewidth=1.2, ls='-',
                            label=r'${\hat{\Psi}}_{MAP}\ =\ %.2f$' %(map_estimate))
	if value_to_label:
	    l = plt.axvline(x=value_to_label, color='r', linewidth=1.2, ls='-', label=label)
            
	# plot credible intervals if given
	if plot_intervals:
	    interval_c1, interval_c2 = ht.compute_credible_intervals(samples_to_plot, plot_intervals)
	    plt.axvline(x=interval_c1, color='#999999', linewidth=0.7, ls='--',
			label=r'%d' %(plot_intervals*100) + '% CI')
	    plt.axvline(x=interval_c2, color='#999999', linewidth=0.7, ls='--')

        # mean of posterior
	if plot_mean:
	    sample_mean = mean(samples_to_plot)
	    plt.axvline(x=sample_mean, color='r', linewidth=0.8, label='Mean')
	if with_legend and plot_intervals:
	    if not bbox_coords:
		lg = plt.legend(handletextpad=0.172, borderpad=0.01, labelspacing=.008,
                                handlelength=1.4, loc='best', numpoints=1)
	    else:
		lg = plt.legend(handletextpad=0.172, borderpad=0.01, labelspacing=.008,
                                handlelength=1.4, loc='best', numpoints=1,
				bbox_to_anchor=bbox_coords)
	    lg.get_frame().set_linewidth(0)
	    for t in lg.get_texts():
		t.set_fontsize(8)
	if title:
	    plt.title(title)
	if plots_dir:
	    if not plot_filename:
		plt.savefig(plots_dir + "sampled_psi_hist_%s" %(plot_id))
	    else:
		plt.savefig(plots_dir + plot_filename + '.pdf')
	return curr_axes
	# Plot joint scores as function of number of samples
	#log_joint_fig = plt.figure(figsize=(7,4.5), dpi=300)
	#skip = 15
	#print "Skip of %d when plotting log joint scores" %(skip)
	#plt.plot(arange(0, len(total_log_scores), skip),
	#	 total_log_scores[arange(0, len(total_log_scores), skip)])
	#print "Total log scores plotted: ", len(total_log_scores)
	#plt.xlabel('Number of iterations (lag not shown)')
	#plt.ylabel('Log joint score')
	#plt.savefig(plots_dir + "log_joint_scores_skip%d_%s" %(skip, plot_id))

def parse_sampler_params(header):
    """
    Parse parameters that were used to produce a set of samples.  Build a Gene instance
    out of the Gene it was used on and return that as well.
    """
    if header[0] == '#':
	# strip header start
	header = header[1:]
    fields = header.split('\t')
    isoforms = []
    params = {}
    exon_lens = {}
    for field in fields:
	key, value = field.split('=')
	if key == 'isoforms':
	    isoform_desc = eval(value)
	elif key == 'exon_lens':
	    exon_lens = dict(eval(value))
	else:
	    params[key] = value
    exons = []
    for e, exon_len in exon_lens.iteritems():
	exons.append(Gene.Exon(0, 0 + exon_len - 1, label=e))
    gene = Gene.Gene(isoform_desc, exons)
    return (params, gene)
	    
def load_samples(samples_filename):
    """
    Load a set of samples from a file and build an associated gene.
    """
    samples_data, h = csv2array(samples_filename, skiprows=1,
                                raw_header=True)
    samples = []
    log_scores = []
    for line in samples_data:
	psi_vals = [float(v) for v in line['sampled_psi'].split(',')]
	samples.append(psi_vals)
	log_scores.append(float(line['log_score']))
    params, gene = parse_sampler_params(h[0])
    return (array(samples), array(log_scores), params, gene)

