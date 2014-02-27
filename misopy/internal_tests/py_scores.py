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


