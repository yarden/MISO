
#include "splicing.h"
#include "splicing_error.h"

#include <math.h>

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI   0.918938533204672741780329736406        
/* log(sqrt(2*pi)) == log(2*pi)/2 */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI    0.398942280401432677939946059934        
/* 1/sqrt(2pi) */
#endif

double splicing_i_dnorm(double x, double mu, double sigma, int plog) {
  if (sigma <= 0) { 
    splicing_error("Invalid `sigma' for normal", __FILE__, __LINE__,
		   SPLICING_EINVAL);
  }

  x = (x - mu) / sigma;

  return (plog ?
	  -(M_LN_SQRT_2PI  +  0.5 * x * x + log(sigma)) :
	  M_1_SQRT_2PI * exp(-0.5 * x * x)  /   sigma);
}

double splicing_dnorm(double x, double mu, double sigma) {
  return splicing_i_dnorm(x, mu, sigma, /*log=*/ 0);
}

double splicing_logdnorm(double x, double mu, double sigma) {
  return splicing_i_dnorm(x, mu, sigma, /*log=*/ 1);
}

