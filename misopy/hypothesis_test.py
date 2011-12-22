##
## Bayesian and frequentist approaches to hypothesis testing for MISO
##
from numpy import *
import os
import scipy
from scipy import stats
from scipy.stats import gaussian_kde
from decimal import Decimal
import misopy
from misopy.samples_utils import *
from misopy.credible_intervals import *

#load_samples, format_credible_intervals, \
#     get_samples_dir_filenames, get_isoforms_from_header

#import matplotlib
#import matplotlib.pyplot as plt
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rcParams['ps.useafm'] = True
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rcParams['font.size'] = 10
# Trying this setting
#plt.rcParams['pdf.fonttype'] = 42

class NullPeakedDensity:
    """
    A density peaked on the null hypothesis
    """
    def __init__(self, dataset):
	self.dataset = dataset

    def evaluate(self, point):
	if point[0] == 0:
	    return inf
	else:
	    return 0

class gaussian_kde_set_covariance(stats.gaussian_kde):
    '''
    from Anne Archibald in mailinglist:
    http://www.nabble.com/Width-of-the-gaussian-in-stats.kde.gaussian_kde---td19558924.html#a19558924
    '''
    def __init__(self, dataset, covariance):
        self.covariance = covariance
        scipy.stats.gaussian_kde.__init__(self, dataset)

    def _compute_covariance(self):
        self.inv_cov = np.linalg.inv(self.covariance)
        self._norm_factor = sqrt(np.linalg.det(2*np.pi*self.covariance)) * self.n
        
class gaussian_kde_covfact(stats.gaussian_kde):
    def __init__(self, dataset, covfact = 'scotts'):
        self.covfact = covfact
        scipy.stats.gaussian_kde.__init__(self, dataset)

    def _compute_covariance_(self):
        '''not used'''
        self.inv_cov = np.linalg.inv(self.covariance)
        self._norm_factor = sqrt(np.linalg.det(2*np.pi*self.covariance)) * self.n

    def covariance_factor(self):
        if self.covfact in ['sc', 'scotts']:
            return self.scotts_factor()
        if self.covfact in ['si', 'silverman']:
            return self.silverman_factor()
        elif self.covfact:
            return float(self.covfact)
        else:
            raise ValueError, \
                'covariance factor has to be scotts, silverman or a number'

    def reset_covfact(self, covfact):
        self.covfact = covfact
        self.covariance_factor()
        self._compute_covariance()
	
def compute_prior_proportion_diff(num_samples):
    """
    Compute the posterior on the difference between two independent proportions (from two
    distinct conditions.)

    The model assumed here is:

    p_c1 ~ Beta(1, 1)
    p_c2 ~ Beta(1, 1)

    We compute P(delta), where delta = p1_c1 - p2_c2, and return a set of num_samples-many samples.
    """
    samples = []
    for n in range(num_samples):
	# sample probabilities from uniform prior
	prob_c1 = random.beta(1, 1)
	prob_c2 = random.beta(1, 1)
	delta = prob_c1 - prob_c2
	samples.append(delta)
    return array(samples)



def compute_delta_densities(samples1_filename, samples2_filename, diff_range,
                            smoothing_param=0.3):
    """
    Compute the Gaussian kernel density fitted distributions over delta for
    the two sets of posterior samples filenames given.  Returns the posterior density
    and prior density as well, assuming a uniform prior over the Psi of the samples
    in the two conditions.
    """
    densities = {}
    # Compute analytic prior density
    prior_density_fn = lambda x: 1 + x if x <= 0 else 1 - x
    analytic_prior_density = map(prior_density_fn, diff_range)
    # Load posterior samples from files
    samples1_results = load_samples(samples1_filename)
    posterior_samples1 = samples1_results[0]
    samples2_results = load_samples(samples2_filename)
    posterior_samples2 = samples2_results[0]

    num_samples, num_isoforms = shape(posterior_samples1)

    # Extract isoforms header information
    header = samples1_results[1]
    isoforms_field = get_isoforms_from_header(header[0])

    # Extract counts from headers
    sample1_counts_info = samples1_results[5]
    sample2_counts_info = samples2_results[5]

    # Record isoform information and counts
    densities['isoforms'] = isoforms_field
    densities['sample1_counts'] = sample1_counts_info
    densities['sample2_counts'] = sample2_counts_info

    # Set prior density function
    densities['prior_density'] = prior_density_fn

    # Posterior samples from MCMC
    densities['samples1'] = posterior_samples1
    densities['samples2'] = posterior_samples2
    
    # Collection of posterior densities (only 1 in two-isoform case)
    densities['posterior_density'] = []

    # Collection of Bayes factors (only 1 in two-isoform case)
    densities['bayes_factor'] = []
    
    # For each isoform, compute its Bayes factor and delta posterior
    warning_outputted = False
    
    for iso_num in range(num_isoforms):
        posterior_diff = posterior_samples1[:, iso_num] - \
                         posterior_samples2[:, iso_num]

        # If the average difference is 0, don't try to fit a kernel to it
        mean_abs_posterior_diff = mean(abs(posterior_diff))

        # If the posterior differences are all identical, the sampler
        # was probably unable to explore the space
        all_same_diff = all(posterior_diff - posterior_diff[0] == 0)

        if all_same_diff and not warning_outputted:
            print "Warning: %s or %s were not properly sampled." \
                  %(samples1_filename, samples2_filename)
            warning_outputted = True
            densities['bayes_factor'] = 0

        if mean_abs_posterior_diff <= .009 or all_same_diff:
            posterior_density = NullPeakedDensity(posterior_diff)
        else:
            # Smoothing by fitting Gaussian kernel density estimator
            posterior_density = gaussian_kde_covfact(posterior_diff, smoothing_param)

        bayes_factor, diff_prior, diff_posterior = compute_bayes_factor(analytic_prior_density,
                                                                        posterior_density)
        densities['bayes_factor'].append(bayes_factor)
        densities['posterior_density'].append(posterior_density)
    
    return densities


def output_samples_comparison(sample1_dir, sample2_dir, output_dir,
                              alpha=.95):
    """
    Compute the bayes factors, posterior means, and other statistics
    between the two samples and output them to a directory.

    Expects two directories with samples from a MISO run, where corresponding
    events in the two samples' directories begin with the same event name.
    """
    print "Given output dir: ", output_dir
    # Retrieve only the files that are in the two given directories
    sample1_filenames = get_samples_dir_filenames(sample1_dir)
    sample2_filenames = get_samples_dir_filenames(sample2_dir)
    
    print "Computing sample comparison between %s and %s..." %(sample1_dir,
                                                               sample2_dir)
    print "  - # files in %s: %d" %(sample1_dir, len(sample1_filenames))
    print "  - # files in %s: %d" %(sample2_dir, len(sample2_filenames))

    # Output header for Bayes factor file
    sample1_label = os.path.basename(os.path.normpath(sample1_dir))
    sample2_label = os.path.basename(os.path.normpath(sample2_dir))
    output_dir = os.path.join(output_dir, "%s_vs_%s" %(sample1_label,
                                                       sample2_label))

    print "Creating comparisons parent directory: %s" %(output_dir)

    # Create parent directory for comparison
    if not os.path.isdir(output_dir):
	os.makedirs(output_dir)
	
    # Create directory for Bayes factors
    bf_output_dir = os.path.join(output_dir, 'bayes-factors/')
    if not os.path.isdir(bf_output_dir):
	os.mkdir(bf_output_dir)
	
    # Create directory for raw delta posteriors
    dp_output_dir = os.path.join(output_dir, 'delta-posteriors/')
    if not os.path.isdir(dp_output_dir):
	os.makedirs(dp_output_dir)
    
    header_fields = ['event_name',
		     'sample1_posterior_mean',
		     'sample1_ci_low',
		     'sample1_ci_high',
		     'sample2_posterior_mean',
		     'sample2_ci_low',
		     'sample2_ci_high',
		     'diff',
		     'bayes_factor',
                     'isoforms',
                     'sample1_counts',
                     'sample1_assigned_counts',
                     'sample2_counts',
                     'sample2_assigned_counts']
    header_line = "\t".join(header_fields) + "\n"
    output_filename = os.path.join(bf_output_dir, "%s_vs_%s.miso_bf" %(sample1_label,
                                                                       sample2_label))
    output_file = open(output_filename, 'w')
    output_file.write(header_line)

    num_events_compared = 0

    # Number of events to put into directories -- used to
    # split up the raw delta-posteriors
    batch_size = 500

    file_num = 0
    curr_batch = file_num

    # Compute the Bayes factors for each file
    for sample1_filename in sample1_filenames:
        split_id = ".miso"
	sample1_event_name = os.path.basename(sample1_filename).split(split_id)[0]
        
	# Find corresponding event filename in sample 2
        sample2_filename = filter(lambda filename:
                                  os.path.basename(filename).split(split_id)[0] == sample1_event_name,
                                  sample2_filenames)

	if len(sample2_filename) == 0:
            continue

        sample2_filename = sample2_filename[0]
        
        num_events_compared += 1
	
	# Compute delta of posterior samples and Bayes factors
	diff_range = arange(-1, 1, 0.001)
	delta_densities = compute_delta_densities(sample1_filename, sample2_filename,
                                                  diff_range)
	bf = delta_densities['bayes_factor']
        
        num_isoforms = shape(delta_densities['samples1'])[1]
        
	sample1_posterior_mean = mean(delta_densities['samples1'], 0)
	sample2_posterior_mean = mean(delta_densities['samples2'], 0)


        # Get the labels of the isoforms
        isoforms_field = delta_densities['isoforms']

        # Get the counts information about both samples
        sample1_counts_info = delta_densities['sample1_counts']
        sample2_counts_info = delta_densities['sample2_counts']

	# Compute posterior mean and credible intervals for sample 1
        sample1_cred_intervals = format_credible_intervals(sample1_event_name,
                                                           delta_densities['samples1'],
                                                           confidence_level=alpha)
        sample1_ci_low = sample1_cred_intervals[2]
        sample1_ci_high = sample1_cred_intervals[3]

        # Compute posterior mean and credible intervals for sample 2
        sample2_cred_intervals = format_credible_intervals(sample1_event_name,
                                                           delta_densities['samples2'],
                                                           confidence_level=alpha)
        sample2_ci_low = sample2_cred_intervals[2]
        sample2_ci_high = sample2_cred_intervals[3]
        
	posterior_diff = sample1_posterior_mean - sample2_posterior_mean

	# Use precision of two decimal places
        if num_isoforms == 2:
            sample1_posterior_mean = Decimal(str(sample1_posterior_mean[0])).quantize(Decimal('0.01'))
            sample2_posterior_mean = Decimal(str(sample2_posterior_mean[0])).quantize(Decimal('0.01'))
            posterior_diff = "%.2f" %(sample1_posterior_mean - sample2_posterior_mean)
            bayes_factor = "%.2f" %(bf[0])
        else:
            posterior_diff = ",".join(["%.2f" %(v) for v in (sample1_posterior_mean - sample2_posterior_mean)])
            sample1_posterior_mean = sample1_cred_intervals[1]
            sample2_posterior_mean = sample2_cred_intervals[1]
            bayes_factor = ",".join(["%.2f" %(max(v, 0)) for v in bf])

        # Write comparison output line
	output_fields = [sample1_event_name,
                         # Mean and confidence bounds for sample 1
			 "%s" %(sample1_posterior_mean),
			 "%s" %(sample1_ci_low),
			 "%s" %(sample1_ci_high),
                         # Mean and confidence bounds for sample 2
			 "%s" %(sample2_posterior_mean),
			 "%s" %(sample2_ci_low),
			 "%s" %(sample2_ci_high),
                         # Delta Psi value
			 "%s" %(posterior_diff),
                         # Bayes factor
			 "%s" %(bayes_factor),
                         # Description of the isoforms
                         "%s" %(isoforms_field),
                         # Counts information for sample 1
                         "%s" %(sample1_counts_info['counts']),
                         "%s" %(sample1_counts_info['assigned_counts']),
                         # Counts information for sample 2
                         "%s" %(sample2_counts_info['counts']),
                         "%s" %(sample2_counts_info['assigned_counts'])]
	output_line = "%s\n" %("\t".join(output_fields))
	output_file.write(output_line)
        
	# Output raw delta posteriors
	dp_header = "delta_posteriors\n"

        # Move to next batch if needed
        if file_num % batch_size == 0:
            curr_batch += 1
            print "Outputting batch number %d (batch size = %d)..." \
                  %(curr_batch, batch_size)

            # Make output dir for the current batch
            batch_dir_name = "batch_%d_%d" %(batch_size, curr_batch)
            curr_dp_dir = os.path.join(dp_output_dir, batch_dir_name)
            
            if not os.path.isdir(curr_dp_dir):
                print "Making output directory: %s" %(curr_dp_dir)
                os.makedirs(curr_dp_dir)
            
        file_num += 1

        # File name for delta posterior file
        dp_filename = os.path.join(curr_dp_dir, sample1_event_name + '.miso_dp')

        # Output the raw delta posteriors
	dp_file = open(dp_filename, 'w')
	dp_file.write(dp_header)

        delta_posteriors = delta_densities['samples1'] - \
                           delta_densities['samples2']

	for delta_posterior in delta_posteriors:
            if num_isoforms == 2:
                delta_posterior = delta_posterior[0:-1]
            dp_output_line = "%s\n" %(",".join(["%.4f" %(v) for v in delta_posterior]))
	    dp_file.write(dp_output_line)
	dp_file.close()
                                 

    print "Compared a total of %d events." %(num_events_compared)
    output_file.close()
    

def compute_bayes_factor(prior_density, posterior_density, at_point=0, print_bayes=False):
    """
    Compute Bayes factor for given fitted densities.
    """
    max_bf = 1e12
    
    # assume prior density is known analytically at delta = 0
    if at_point == 0:
	diff_prior = 1
    else:
	diff_prior = prior_density.evaluate([at_point])
    diff_posterior = posterior_density.evaluate([at_point])

    if diff_posterior == 0:
	bayes_factor = max_bf
    elif diff_posterior == inf:
	bayes_factor = 0
    else:
	# Compute factor relative to alternative hypothesis
	bayes_factor = diff_prior / diff_posterior
        bayes_factor = bayes_factor[0]

    if print_bayes:
	print "diff_posterior: %.4f" %(diff_posterior)
	print "bayes_factor: %.2f" %(bayes_factor)
	
    return bayes_factor, diff_prior, diff_posterior
    
def main():
    pass
    
if __name__ == '__main__':
    main()

