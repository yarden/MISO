##
## Pure Python scoring functions as reference implementation
##
## Yarden Katz <yarden@mit.edu>
##
import os
import sys
import time

import numpy as np
import numpy.linalg as linalg
from numpy.random import multivariate_normal
from numpy.random import normal

import math

import scipy
import scipy.misc 
from scipy.special import gammaln

def set_diag(a, v):
    for i, elt in enumerate(a):
        a[i, i] = v
    return a


def dirichlet_lnpdf(alpha, x):
    """
    Substitute for dirichlet_lnpdf of pygsl.
    """
    dir_log_pdf = \
        gammaln(np.sum(alpha)) - sum(gammaln(alpha)) + np.dot((alpha - 1).T, np.log(x).T)
    return dir_log_pdf


def approx_eq(p1, p2, error=0.0001):
    return (np.abs(p1 - p2) < error)


def original_logistic_normal_log_pdf(theta, mu, sigma):
    """
    The log of the PDF for the multivariate Logistic-Normal distribution.
    Written in terms of k-1 dimensions.
    """
    theta = theta[:-1]
    if len(theta) != len(mu):
        raise Exception, "len(theta) = %d != len(mu) = %d" \
              %(len(theta), len(mu))
    theta_last = 1-sum(theta)
    covar_constant = np.power(float(linalg.det(2*math.pi*sigma)), -0.5)
    invsigma = linalg.inv(np.transpose(sigma))
    prod_theta = 1/(float(np.prod(theta))*theta_last)
    first_log = -0.5*np.transpose(np.log(theta/theta_last) - mu)
    second_log = np.log(theta/theta_last) - mu
    exp_part = np.dot(np.dot(first_log, invsigma), second_log)
    pdf_value = covar_constant*prod_theta*np.exp(exp_part)
    return np.log(pdf_value)


def propose_norm_drift_psi_alpha(alpha_vector,
                                 sigma_proposal):
    if len(alpha_vector) == 1:
        alpha_vector = alpha_vector[0]
#            print "proposing from normal with mean: ", alpha_vector, " exp: ", exp(alpha_vector)
        alpha_next = [normal(alpha_vector, sigma_proposal)]
#            print "got alpha_next: ", alpha_next, " exp: ", exp(alpha_next)
        new_psi = logit_inv([alpha_next[0]])[0]
        new_psi_vector = [new_psi, 1-new_psi]
    else:
        alpha_next = multivariate_normal(alpha_vector, sigma_proposal)
        new_psi = logit_inv(alpha_next)
        new_psi_vector = np.concatenate((new_psi, array([1-sum(new_psi)])))
    return (new_psi_vector, alpha_next)


def logit_inv(x):
    """
    Takes a value on x \in (-inf, inf) and transforms it to a value on (0, 1).
    """
    #p = exp(x)/(1+exp(x))
    denom = np.lib.function_base.append(x, 0)
    p = np.exp(x)/(np.sum(np.exp(denom)))
    return p


def sample_from_multivar_normal(mu, cov):
    # Define mean, covariance matrix, and no. of samples required
    Ndraws = 1000

    # Do factorisation (can store this and use it again later)
    L = np.linalg.cholesky(cov)

    # Get 3*Ndraws Gaussian random variables (mean=0, variance=1)
    norm = np.random.normal(size=Ndraws*3).reshape(3, Ndraws)

    # Construct final set of random numbers (with correct mean)
    rand = mean + np.dot(L, norm)




