from pysplicing import *
MISO_START_AUTO=0L
MISO_START_UNIFORM=1L
MISO_START_RANDOM=2L
MISO_START_GIVEN=3L
MISO_START_LINEAR=4L
# stop after fixed no. iterations
MISO_STOP_FIXEDNO=0L
MISO_STOP_CONVERGENT_MEAN=1L
# various algorithms
MISO_ALGO_REASSIGN=0L
MISO_ALGO_MARGINAL=1L
MISO_ALGO_CLASSES=2L

# priors
MISO_PRIOR_AUTO=0L
MISO_PRIOR_DIRICHLET=1L
MISO_PRIOR_LOGISTIC=2L

def doMISO(GFF, gene, read_pos, read_cigar, read_len, num_samples,
           burn_in, lag, prior = MISO_PRIOR_DIRICHLET,
           dirichlet_prior_params = None, logistic_prior_mean = 0.0,
           logistic_prior_var = 3.0):

    if dirichlet_prior_params is None:
        no_iso = noIso(GFF)[0]
        dirichlet_prior_params = (1,) * no_iso

    return pysplicing.MISO(
        GFF, gene, read_pos, read_cigar, read_len, num_samples, burn_in,
        lag, dirichlet_prior_params)
