cimport cython

from cython.view cimport array as cvarray
from cpython.array cimport array, clone

cimport matrix_utils

from libc.math cimport log
from libc.math cimport exp
from libc.math cimport pow
from libc.stdlib cimport rand

# Math constants
from libc.math cimport M_PI

cdef float MY_MAX_INT = float(10000)


cpdef double[:] \
  my_cumsum(double[:] input_array,
            double[:] cumsum_array):
    """
    Return cumulative sum of array.
    """
    cdef int num_elts = input_array.shape[0]
    # Cumulative sum at every position
    cdef int curr_elt = 0
    # Current cumulative sum: starts at first element
    cdef double curr_cumsum = 0.0
    for curr_elt in xrange(num_elts):
        cumsum_array[curr_elt] = (input_array[curr_elt] + curr_cumsum)
        curr_cumsum = cumsum_array[curr_elt]
    return cumsum_array


cpdef double dirichlet_lnpdf(double[:] alpha,
                             double[:] vector):
    """
    Wrapper for dirichlet log pdf scoring function.
    """
    cdef int D = vector.size
    return dirichlet_log_pdf_raw(D,
                                 &alpha[0], alpha.strides[0],
                                 &vector[0], vector.strides[0])

# from Borg project
@cython.infer_types(True)
cdef double dirichlet_log_pdf_raw(int D,
                                  double* alpha,
                                  int alpha_stride,
                                  double* vector,
                                  int vector_stride,):
    """Compute the log of the Dirichlet PDF evaluated at one vector."""
    cdef void* alpha_p = alpha
    cdef void* vector_p = vector

    # first term
    term_a = 0.0

    for d in xrange(D):
        term_a += (<double*>(alpha_p + alpha_stride * d))[0]

    term_a = libc.math.lgamma(term_a)

    # second term
    term_b = 0.0

    for d in xrange(D):
        term_b += libc.math.lgamma((<double*>(alpha_p + alpha_stride * d))[0])

    # third term
    cdef double alpha_d
    cdef double vector_d

    term_c = 0.0

    for d in xrange(D):
        alpha_d = (<double*>(alpha_p + alpha_stride * d))[0]
        vector_d = (<double*>(vector_p + vector_stride * d))[0]

        term_c += (alpha_d - 1.0) * libc.math.log(vector_d)

    # ...
    return term_a - term_b + term_c



# def logistic_normal_log_pdf(self, theta, mu):
#     """
#     The log of the PDF for the multivariate Logistic-Normal distribution.
#     Written in terms of k-1 dimensions.
#     """
#     theta = theta[:-1]
#     if len(theta) != len(mu):
#         raise Exception, "len(theta) = %d != len(mu) = %d" \
#               %(len(theta), len(mu))
#     theta_last = 1-sum(theta)
#     sigma = self.params['sigma_proposal']
#     covar_constant = power(float(la.det(2*math.pi*sigma)), -0.5)
#     invsigma = la.inv(transpose(sigma))
#     prod_theta = 1/(float(prod(theta))*theta_last)
#     first_log = -0.5*transpose(log(theta/theta_last) - mu)
#     second_log = log(theta/theta_last) - mu
#     exp_part = dot(dot(first_log, invsigma), second_log)
#     pdf_value = covar_constant*prod_theta*exp(exp_part)
#     return log(pdf_value)
        

# cdef double logistic_normal_log_pdf(np.ndarray[double, ndim=1] theta,
#                                     np.ndarray[double, ndim=1] mu,
#                                     np.ndarray[double, ndim=2] sigma):
#     """
#     The log of the PDF for the multivariate Logistic-Normal distribution.
#     Written in terms of k-1 dimensions.

#     theta : k dimensional vector to score
#     mu : mean parameter
#     sigma : covariance matrix (2d array)
#     k : dimension
#     """
#     cdef double result = 0.0
#     # Use k-1 dimensional theta vector
#     theta = theta[:-1]
#     assert theta.shape[0] == mu.shape[0], \
#       "Length of theta and mu do not match."
#     theta_last = 1 - matrix_utils.sum_array(theta, theta.shape[0])
#     covar_constant = np.power(float(linalg.det(2*math.pi*sigma)), -0.5)
#     invsigma = linalg.inv(np.transpose(sigma))
#     prod_theta = 1/(float(np.prod(theta))*theta_last)
#     first_log = -0.5*np.transpose(log(theta/theta_last) - mu)
#     second_log = log(theta/theta_last) - mu
#     exp_part = np.dot(np.dot(first_log, invsigma), second_log)
#     pdf_value = covar_constant*prod_theta*np.exp(exp_part)
#     return log(pdf_value)

cpdef double logistic_normal_log_pdf(double[:] theta,
                                     double[:] mu,
                                     double sigma):
    """
    The log of the PDF for the multivariate Logistic-Normal distribution.
    Assumes that Sigma is a diagonal matrix where the diagonal
    is all set to the constant 'sigma'

    theta : k dimensional vector to score
    mu : mean parameter
    sigma : float value that is the diagonal of Sigma covar matrix
    k : dimension
    """
    cdef double ltheta = 1.0
    cdef double prodTheta = 1.0
    cdef double expPart = 0.0
    cdef double pdfVal
    cdef int vect_len = theta.shape[0]
    cdef double covarConst = pow(2 * M_PI * sigma, -0.5 * vect_len)
    cdef int i
    cdef double at
    cdef double tmp

    for i in xrange(vect_len):
        at = theta[i]
        ltheta -= at
        prodTheta *= at

    prodTheta = 1.0 / prodTheta / ltheta

    for i in xrange(vect_len):
        tmp = log(theta[i] / ltheta) - mu[i]
        expPart += (-0.5) * tmp * tmp / sigma
    pdfVal = covarConst * prodTheta * exp(expPart)
    score = log(pdfVal)
    return score


# int splicing_mvplogisnorm(const splicing_vector_t *theta, 
# 			  const splicing_vector_t *mu, 
# 			  double sigma, int len, double *score) {
  
#   double covarConst = pow(2 * M_PI * sigma, -0.5 * len);
#   double ltheta=1.0, prodTheta=1.0, expPart=0.0, pdfVal;
#   int i;
  
#   for (i=0; i<len; i++) {
#     double at=VECTOR(*theta)[i];
#     ltheta -= at;
#     prodTheta *= at;
#   }
#   prodTheta = 1.0 / prodTheta / ltheta;
  
#   for (i=0; i<len; i++) {
#     double tmp=log(VECTOR(*theta)[i] / ltheta) - VECTOR(*mu)[i];
#     expPart += (-0.5) * tmp * tmp / sigma;
#   }
  
#   pdfVal = covarConst * prodTheta * exp(expPart);

#   *score = log(pdfVal);
  
#   return 0;
# }


# from libc.math cimport pow

# @cython.boundscheck(False)
# @cython.wraparound(False)
# cpdef double determinant(double[:, :] M):
#     '''
#     Compute the determinant of a square nxn matrix using
#     the recursive formula:
 
#         det(M) = sum_{i=0}^{n-1} (-1^i)*M[0,i]*det(M_minor(0,i))

#     where M_minor(0,i) is the (n-1)x(n-1) matrix formed by
#     removing row 0 and column i from M.
    
#     source: http://nbviewer.ipython.org/github/carljv/cython_testing/blob/master/cython_linalg.ipynb
#     '''
#     assert M.shape[0] == M.shape[1], 'Matrix is not square.'
    
#     cdef int i, j
#     cdef int n = M.shape[0]
#     cdef double det = 0.0
#     cdef double coef
#     cdef double[:, :] M_minor = np.empty((n-1, n-1))
    
#     if n == 1:
#         # If M is a scalar (1x1) just return it
#         return M[0, 0]
#     else:
#         # If M is nxn, then get its (n-1)x(n-1) minors
#         # (one for each of M's n columns) and compute 
#         # their determinants and add them to the summation.
#         for j in xrange(n):a
#             coef = pow(-1, j) * M[0, j]
#             _get_minor(M, M_minor, 0, j)
#             det += coef * determinant(M_minor)
#         return det

    
# @cython.boundscheck(False)
# @cython.wraparound(False)
# cdef void _get_minor(double[:, :] M, double[:, :] M_minor, 
#                     int row, int col):
#     '''
#     Return the minor of a matrix, by removing a specified
#     row and column

#     If M is nxn, then _get_minor(M, row, col) will fill in
#     the (n-1)x(n-1) matrix M_minor by removing row `row` 
#     and column `col` from M.
    
#     source: http://nbviewer.ipython.org/github/carljv/cython_testing/blob/master/cython_linalg.ipynb
#     '''
#     cdef:
#         int n = M.shape[0]
#         int i_to, j_to, i_from, j_from
  
        
#     # _from indicates the index of the original
#     # matrix M, _to, indicates the index of the
#     # result matrix M_minor.
#     i_from = 0
#     for i_to in xrange(n-1):
#         if (i_to == row): 
#             # This is the row to exclude from the
#             # minor, so skip it.
#             i_from += 1
#         j_from = 0
#         for j_to in xrange(n-1):
#             if (j_to == col): 
#             # This is the column to exclude from the
#             # minor, so skip it.
#                 j_from += 1
            
#             M_minor[i_to, j_to] = M[i_from, j_from]
#             j_from += 1
        
#         i_from += 1
