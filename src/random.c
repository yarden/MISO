/* -*- mode: C -*-  */

#include "splicing_random.h"
#include "splicing_error.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <limits.h>
#include <string.h>
#include "splicing_vector.h"
#include "splicing_memory.h"
#include "splicing.h"

/** 
 * \section about_rngs 
 * 
 * <section>
 * <title>About random numbers in splicing, use cases</title>
 * 
 * <para> 
 * Some algorithms in splicing, e.g. the generation of random graphs,
 * require random number generators (RNGs). Prior to version 0.6
 * splicing did not have a sophisticated way to deal with random number
 * generators at the C level, but this has changed. From version 0.6
 * different and multiple random number generators are supported.
 * </para>
 * </section>
 * 
 */

/** 
 * \section rng_use_cases
 * 
 * <section><title>Use cases</title>
 * 
 * <section><title>Normal (default) use</title>
 * <para> 
 * If the user does not use any of the RNG functions explicitly, but calls
 * some of the randomized splicing functions, then a default RNG is set
 * up the first time an splicing function needs random numbers. The
 * seed of this RNG is the output of the <code>time(0)</code> function
 * call, using the <code>time</code> function from the standard C
 * library. This ensures that splicing creates a different random graph,
 * each time the C program is called.
 * </para>
 * 
 * <para> 
 * The created default generator is stored in the \ref
 * splicing_rng_default variable.
 * </para>
 * </section>
 * 
 * <section><title>Reproducible simulations</title>
 * <para> 
 * If reproducible results are needed, then the user should set the
 * seed of the default random number generator explicitly, using the 
 * \ref splicing_rng_seed() function on the default generator, \ref
 * splicing_rng_default. When setting the seed to the same number,
 * splicing generates exactly the same random graph (or series of random
 * graphs).
 * </para>
 * </section>
 * 
 * <section><title>Changing the default generator</title>
 * <para> 
 * By default splicing uses the \ref splicing_rng_default random number
 * generator. This can be changed any time by calling \ref
 * splicing_rng_set_default(), with an already initialized random number
 * generator. Note that the old (replaced) generator is not
 * destroyed, so no memory is deallocated.
 * </para>
 * </section>
 *
 * <section><title>Using multiple generators</title>
 * <para> 
 * splicing also provides functions to set up multiple random number
 * generators, using the \ref splicing_rng_init() function, and then
 * generating random numbers from them, e.g. with \ref splicing_rng_get_integer()
 * and/or \ref splicing_rng_get_unif() calls. 
 * </para>
 * 
 * <para>
 * Note that initializing a new random number generator is
 * independent of the generator that the splicing functions themselves
 * use. If you want to replace that, then please use \ref
 * splicing_rng_set_default().
 * </para>
 * </section>
 *
 * <section><title>Example</title>
 * <para>
 * \example examples/simple/random_seed.c
 * </para>
 * </section>
 *
 * </section>
 */

/* ------------------------------------ */

/* This is a placeholder, essentially. The R RNG is used from R and 
   the Python RNG is used from Python. */

splicing_rng_type_t splicing_rngtype_null = { 
  /* name= */      "null",
  /* min=  */      0,
  /* max=  */      0,
  /* init= */      0,
  /* destroy= */   0,
  /* seed= */      0,
  /* get= */       0,
  /* get_real= */  0,
  /* get_norm= */  0,
  /* get_gamma= */ 0
};

typedef struct {
} splicing_i_rng_null_state_t;

/** 
 * \function splicing_rng_set_default
 * Set the default splicing random number generator
 * 
 * \param rng The random number generator to use as default from now
 *    on. Calling \ref splicing_rng_destroy() on it, while it is still
 *    being used as the default will result in crashes and/or
 *    unpredictable results.
 * 
 * Time complexity: O(1).
 */

void splicing_rng_set_default(splicing_rng_t *rng) {
  splicing_rng_default = (*rng);
}

splicing_i_rng_null_state_t splicing_i_rng_default_state;

#ifndef USING_R

#define addr(a) (&a)

/** 
 * \var splicing_rng_default
 * The default splicing random number generator
 * 
 * This generator is used by all builtin splicing functions that need to
 * generate random numbers; e.g. all random graph generators. 
 * 
 * You can use \ref splicing_rng_default with \ref splicing_rng_seed()
 * to set its seed.
 * 
 * You can change the default generator using the \ref
 * splicing_rng_set_default() function. 
 */

splicing_rng_t splicing_rng_default = { 
  addr(splicing_rngtype_null),
  addr(splicing_i_rng_default_state),
  /* def= */ 1
};

#undef addr

#endif

/* ------------------------------------ */

#ifdef USING_R

double  unif_rand(void);
double  norm_rand(void);
double  Rf_rgamma(double, double);


int splicing_rng_R_init(void **state) {
  SPLICING_ERROR("R RNG error, unsupported function called",
	       SPLICING_EINTERNAL);
  return 0;
}

void splicing_rng_R_destroy(void *state) {
  splicing_error("R RNG error, unsupported function called",
	       __FILE__, __LINE__, SPLICING_EINTERNAL);
}

int splicing_rng_R_seed(void *state, unsigned long int seed) {
  SPLICING_ERROR("R RNG error, unsupported function called",
	       SPLICING_EINTERNAL);
  return 0;
}

unsigned long int splicing_rng_R_get(void *state) {
  return unif_rand() * 0x7FFFFFFFUL;
}

double splicing_rng_R_get_real(void *state) {
  return unif_rand();
}

double splicing_rng_R_get_norm(void *state) {
  return norm_rand();
}

double splicing_rng_R_get_gamma(void *state, double a) {
  return Rf_rgamma(a, /*scale=*/ 1.0);
}

splicing_rng_type_t splicing_rngtype_R = {
  /* name= */      "GNU R",
  /* min=  */      0,
  /* max=  */      0x7FFFFFFFUL,
  /* init= */      splicing_rng_R_init,
  /* destroy= */   splicing_rng_R_destroy,
  /* seed= */      splicing_rng_R_seed,
  /* get= */       splicing_rng_R_get,
  /* get_real= */  splicing_rng_R_get_real,
  /* get_norm= */  splicing_rng_R_get_norm,
  /* get_gamma= */ splicing_rng_R_get_gamma
};

splicing_rng_t splicing_rng_default = { 
  &splicing_rngtype_R,
  0,
  /* def= */ 1
};

#endif

/* ------------------------------------ */

double splicing_rgamma(splicing_rng_t *rng, double a);

/** 
 * \function splicing_rng_init
 * Initialize a random number generator
 * 
 * This function allocates memory for a random number generator, with
 * the given type, and sets its seed to the default.
 * 
 * \param rng Pointer to an uninitialized RNG.
 * \param type The type of the RNG, please see the documentation for
 *    the supported types.
 * \return Error code.
 * 
 * Time complexity: depends on the type of the generator, but usually
 * it should be O(1).
 */

int splicing_rng_init(splicing_rng_t *rng, const splicing_rng_type_t *type) {
  rng->type=type;
  SPLICING_CHECK(rng->type->init(&rng->state));
  return 0;
}

/** 
 * \function splicing_rng_destroy
 * Deallocate memory associated with a random number generator
 * 
 * \param rng The RNG to destroy. Do not destroy an RNG that is used
 *    as the default splicing RNG. 
 * 
 * Time complexity: O(1).
 */

void splicing_rng_destroy(splicing_rng_t *rng) {
  rng->type->destroy(rng->state);
}

/**
 * \function splicing_rng_seed
 * Set the seed of a random number generator
 * 
 * \param rng The RNG. 
 * \param seed The new seed.
 * \return Error code.
 * 
 * Time complexity: usually O(1), but may depend on the type of the
 * RNG.
 */
int splicing_rng_seed(splicing_rng_t *rng, unsigned long int seed) {
  const splicing_rng_type_t *type=rng->type;
  rng->def=0;
  SPLICING_CHECK(type->seed(rng->state, seed));
  return 0;
}

/** 
 * \function splicing_rng_max 
 * Query the maximum possible integer for a random number generator
 * 
 * \param rng The RNG.
 * \return The largest possible integer that can be generated by
 *         calling \ref splicing_rng_get_integer() on the RNG.
 * 
 * Time complexity: O(1).
 */

unsigned long int splicing_rng_max(splicing_rng_t *rng) {
  const splicing_rng_type_t *type=rng->type;
  return type->max;
}

/**
 * \function splicing_rng_min
 * Query the minimum possible integer for a random number generator
 * 
 * \param rng The RNG.
 * \return The smallest possible integer that can be generated by
 *         calling \ref splicing_rng_get_integer() on the RNG.
 * 
 * Time complexity: O(1).
 */

unsigned long int splicing_rng_min(splicing_rng_t *rng) {
  const splicing_rng_type_t *type=rng->type;
  return type->min;
}

/** 
 * \function splicing_rng_name
 * Query the type of a random number generator
 * 
 * \param rng The RNG.
 * \return The name of the type of the generator. Do not deallocate or
 *         change the returned string pointer.
 * 
 * Time complexity: O(1).
 */

const char *splicing_rng_name(splicing_rng_t *rng) {
  const splicing_rng_type_t *type=rng->type;
  return type->name;
}

/** 
 * \function splicing_rng_get_integer
 * Generate an integer random number from an interval
 * 
 * \param rng Pointer to the RNG to use for the generation. Use \ref
 *        splicing_rng_default here to use the default splicing RNG.
 * \param l Lower limit, inclusive, it can be negative as well.
 * \param h Upper limit, inclusive, it can be negative as well, but it
 *        should be at least <code>l</code>.
 * \return The generated random integer.
 * 
 * Time complexity: depends on the generator, but should be usually
 * O(1).
 */

long int splicing_rng_get_integer(splicing_rng_t *rng,
				long int l, long int h) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_real) {
    return (long int)(type->get_real(rng->state)*(h-l+1)+l);
  } else if (type->get) {
    unsigned long int max=type->max;
    return (long int)(type->get(rng->state))/
      ((double)max+1)*(h-l+1)+l;
  }
  SPLICING_ERROR("Internal random generator error", SPLICING_EINTERNAL);
  return 0;
}

/** 
 * \function splicing_rng_get_normal
 * Normally distributed random numbers
 * 
 * \param rng Pointer to the RNG to use. Use \ref splicing_rng_default
 *        here to use the default splicing RNG.
 * \param m The mean.
 * \param s Standard deviation.
 * \return The generated normally distributed random number.
 * 
 * Time complexity: depends on the type of the RNG.
 */

double splicing_rng_get_normal(splicing_rng_t *rng, 
				    double m, double s) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_norm) {
    return type->get_norm(rng->state)*s+m;
  } else{
    SPLICING_ERROR("Internal RNG error", SPLICING_EINTERNAL);
  }
}

/** 
 * \function splicing_rng_get_unif
 * Generate real, uniform random numbers from an interval
 * 
 * \param rng Pointer to the RNG to use. Use \ref splicing_rng_default
 *        here to use the default splicing RNG.
 * \param l The lower bound, it can be negative.
 * \param h The upper bound, it can be negative, but it has to be
 *        larger than the lower bound.
 * \return The generated uniformly distributed random number.
 * 
 * Time complexity: depends on the type of the RNG.
 */

double splicing_rng_get_unif(splicing_rng_t *rng, 
				  double l, double h) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_real) {
    return type->get_real(rng->state)*(h-l)+l;
  } else if (type->get) {
    unsigned long int max=type->max;
    return type->get(rng->state)/((double)max+1)*(double)(h-l)+l;
  }
  SPLICING_ERROR("Internal random generator error", SPLICING_EINTERNAL);
  return 0;  
}

/** 
 * \function splicing_rng_get_unif01
 * Generate real, uniform random number from the unit interval
 *
 * \param rng Pointer to the RNG to use. Use \ref splicing_rng_default
 *        here to use the default splicing RNG.
 * \return The generated uniformly distributed random number.
 * 
 * Time complexity: depends on the type of the RNG.
 */

double splicing_rng_get_unif01(splicing_rng_t *rng) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_real) {
    return type->get_real(rng->state);
  } else if (type->get) {
    unsigned long int max=type->max;
    return type->get(rng->state)/((double)max+1);
  }
  SPLICING_ERROR("Internal random generator error", SPLICING_EINTERNAL);
  return 0;  
}

double splicing_rng_get_gamma(splicing_rng_t *rng, double a) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_gamma) {
    return type->get_gamma(rng->state, a);
  } else {
    return splicing_rgamma(rng, a);
  }
}

unsigned long int splicing_rng_get_int31(splicing_rng_t *rng) {
  const splicing_rng_type_t *type=rng->type;
  unsigned long int max=type->max;
  if (type->get && max==0x7FFFFFFFUL) {
    return type->get(rng->state);
  } else if (type->get_real) {
    return type->get_real(rng->state)*0x7FFFFFFFUL;
  } else { 
    return splicing_rng_get_unif01(rng)*0x7FFFFFFFUL;
  }
}

int splicing_rng_inited = 0;

#ifndef HAVE_EXPM1
#ifndef USING_R			/* R provides a replacement */
/* expm1 replacement */
static double expm1 (double x)
{
    if (fabs(x) < M_LN2)
    {
        /* Compute the Taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */

        double i = 1.0;
        double sum = x;
        double term = x / 1.0;

        do
        {
            term *= x / ++i;
            sum += term;
        }
        while (fabs(term) > fabs(sum) * 2.22e-16);
      
        return sum;
    }

    return expl(x) - 1.0L;
}
#endif
#endif

#ifndef HAVE_RINT
#ifndef USING_R			/* R provides a replacement */
/* rint replacement */
static double rint (double x)
{
   return ( (x<0.) ? -floor(-x+.5) : floor(x+.5) );
}
#endif
#endif

#ifndef HAVE_RINTF
static float rintf (float x)
{
   return ( (x<(float)0.) ? -(float)floor(-x+.5) : (float)floor(x+.5) );
}
#endif

/*
 * \ingroup internal
 * 
 * This function appends the rest of the needed random number to the 
 * result vector.
 */

int splicing_random_sample_alga(splicing_vector_t *res, int l, int h, 
			      int length) {
  double N=h-l+1;
  double n=length;
  
  double top=N-n;
  double Nreal=N;
  double S=0;
  double V, quot;
  
  l=l-1;

  while (n>=2) {
    V=RNG_UNIF01();
    S=1;
    quot=top/Nreal;
    while (quot>V) {
      S+=1;
      top=-1.0+top;
      Nreal=-1.0+Nreal;
      quot=(quot*top)/Nreal;
    }
    l+=S;
    splicing_vector_push_back(res, l);	/* allocated */
    Nreal=-1.0+Nreal; n=-1+n;
  }
  
  S=floor(round(Nreal)*RNG_UNIF01());
  l+=S+1;
  splicing_vector_push_back(res, l);	/* allocated */
  
  return 0;
}

/**
 * \ingroup nongraph
 * \function splicing_random_sample
 * \brief Generates an increasing random sequence of integers.
 * 
 * </para><para>
 * This function generates an increasing sequence of random integer
 * numbers from a given interval. The algorithm is taken literally
 * from (Vitter 1987). This method can be used for generating numbers from a
 * \em very large interval. It is primarily created for randomly
 * selecting some edges from the sometimes huge set of possible edges
 * in a large graph.
 * </para><para>
 * Note that the type of the lower and the upper limit is \c double,
 * not \c int. This does not mean that you can pass fractional
 * numbers there; these values must still be integral, but we need the
 * longer range of \c double in several places in the library
 * (for instance, when generating Erdos-Renyi graphs).
 * \param res Pointer to an initialized vector. This will hold the
 *        result. It will be resized to the proper size.
 * \param l The lower limit of the generation interval (inclusive). This must
 *        be less than or equal to the upper limit, and it must be integral.
 *        Passing a fractional number here results in undefined behaviour.
 * \param h The upper limit of the generation interval (inclusive). This must
 *        be greater than or equal to the lower limit, and it must be integral.
 *        Passing a fractional number here results in undefined behaviour.
 * \param length The number of random integers to generate.
 * \return The error code \c SPLICING_EINVAL is returned in each of the
 *         following cases: (1) The given lower limit is greater than the
 *         given upper limit, i.e. \c l &gt; \c h. (2) Assuming that
 *         \c l &lt; \c h and N is the sample size, the above error code is
 *         returned if N &gt; |\c h - \c l|, i.e. the sample size exceeds the
 *         size of the candidate pool.
 *
 * Time complexity: according to (Vitter 1987), the expected
 * running time is O(length).
 *
 * </para><para>
 * Reference:
 * \clist
 * \cli (Vitter 1987)
 *   J. S. Vitter. An efficient algorithm for sequential random sampling.
 *   \emb ACM Transactions on Mathematical Software, \eme 13(1):58--67, 1987.
 * \endclist
 *
 * \example examples/simple/splicing_random_sample.c
 */

int splicing_random_sample(splicing_vector_t *res, double l, double h, 
			 int length) {
  double N=h-l+1;
  double n=length;
  int retval;

  double nreal=length;
  double ninv=1.0/nreal;
  double Nreal=N;
  double Vprime;
  double qu1=-n+1+N;
  double qu1real=-nreal+1.0+Nreal;
  double negalphainv=-13;
  double threshold=-negalphainv*n;
  double S;

  /* getting back some sense of sanity */
  if (l > h)
    SPLICING_ERROR("Lower limit is greater than upper limit", SPLICING_EINVAL);
  /* now we know that l <= h */
  if (length > N)
    SPLICING_ERROR("Sample size exceeds size of candidate pool", SPLICING_EINVAL);

  /* treat rare cases quickly */
  if (l==h) {
    SPLICING_CHECK(splicing_vector_resize(res, 1));
    VECTOR(*res)[0] = l;
    return 0;
  }
  if (length==N) {
    long int i = 0;
    SPLICING_CHECK(splicing_vector_resize(res, length));
    for (i = 0; i < length; i++) {
      VECTOR(*res)[i] = l++;
    }
    return 0;
  }

  splicing_vector_clear(res);
  SPLICING_CHECK(splicing_vector_reserve(res, length));

  Vprime=exp(log(RNG_UNIF01())*ninv);
  l=l-1;

  while (n>1 && threshold < N) {
    double X, U;
    double limit, t;
    double negSreal, y1, y2, top, bottom;
    double nmin1inv=1.0/(-1.0+nreal);
    while (1) {
      while(1) {
	X=Nreal*(-Vprime+1.0);
	S=floor(X);
	// if (S==0) { S=1; }
	if (S <qu1) { break; }
	Vprime = exp(log(RNG_UNIF01())*ninv);
      }
      U=RNG_UNIF01();
      negSreal=-S;
      
      y1=exp(log(U*Nreal/qu1real)*nmin1inv);
      Vprime=y1*(-X/Nreal+1.0)*(qu1real/(negSreal+qu1real));
      if (Vprime <= 1.0) { break; }
      
      y2=1.0;
      top=-1.0+Nreal;
      if (-1+n > S) {
	bottom=-nreal+Nreal; 
	limit=-S+N;
      } else {
	bottom=-1.0+negSreal+Nreal;
	limit=qu1;
      }
      for (t=-1+N; t>=limit; t--) {
	y2=(y2*top)/bottom;
	top=-1.0+top;
	bottom=-1.0+bottom;
      }
      if (Nreal/(-X+Nreal) >= y1*exp(log(y2)*nmin1inv)) {
	Vprime=exp(log(RNG_UNIF01())*nmin1inv);
	break;
      }
      Vprime=exp(log(RNG_UNIF01())*ninv);
    }
        
    l+=S+1;
    splicing_vector_push_back(res, l);	/* allocated */
    N=-S+(-1+N);   Nreal=negSreal+(-1.0+Nreal);
    n=-1+n;   nreal=-1.0+nreal; ninv=nmin1inv;
    qu1=-S+qu1; qu1real=negSreal+qu1real;
    threshold=threshold+negalphainv;
  }
  
  if (n>1) {
    retval=splicing_random_sample_alga(res, l+1, h, n);
  } else {
    retval=0;
    S=floor(N*Vprime);
    l+=S+1;
    splicing_vector_push_back(res, l);	/* allocated */
  }

  return retval;
}

/**
 * \ingroup nongraph
 * \function splicing_fisher_yates_shuffle
 * \brief The Fisher-Yates shuffle, otherwise known as Knuth shuffle.
 *
 * The Fisher-Yates shuffle generates a random permutation of a sequence.
 *
 * \param seq A sequence of at least one object, indexed from zero. A
 *        sequence of one object is trivially permuted, hence nothing further
 *        is done with this sequence.
 * \return Error code:
 *         \clist
 *         \cli SPLICING_EINVAL
 *           Invalid sequence; \p seq is a null object or it has zero elements.
 *         \endclist
 *
 * Time complexity: depends on the random number generator used, but should
 * usually be O(n), where n is the length of \p seq.
 *
 * </para><para>
 * References:
 * \clist
 * \cli (Fisher &amp; Yates 1963)
 *   R. A. Fisher and F. Yates. \emb Statistical Tables for Biological,
 *   Agricultural and Medical Research. \eme Oliver and Boyd, 6th edition,
 *   1963, page 37.
 * \cli (Knuth 1998)
 *   D. E. Knuth. \emb Seminumerical Algorithms, \eme volume 2 of \emb The Art
 *   of Computer Programming. \eme Addison-Wesley, 3rd edition, 1998, page 145.
 * \endclist
 *
 * \example examples/simple/splicing_fisher_yates_shuffle.c
 */

int splicing_fisher_yates_shuffle(splicing_vector_t *seq) {
  /* sanity checks */
  if (seq == NULL) {
    SPLICING_ERROR("Sequence is a null pointer", SPLICING_EINVAL);
  }
  if (splicing_vector_size(seq) < 1) {
    SPLICING_ERROR("Empty sequence", SPLICING_EINVAL);
  }
  /* obtain random permutation */
  long int i, k;
  for (i = splicing_vector_size(seq) - 1; i > 0; i--) {
    k = RNG_INTEGER(0, i);  /* 0 <= k <= i */
    /* We possibly have k == i, in which case we leave seq[i] alone. */
    SPLICING_CHECK(splicing_vector_swap_elements(seq, k, i));
  }

  return 0;
}
  
#ifdef USING_R

/* These are never called. */

double splicing_norm_rand(splicing_rng_t *rng) {
  return norm_rand();
}

double splicing_rgeom(splicing_rng_t *rng, double p) {
  return Rf_rgeom(p);
}

double splicing_rbinom(splicing_rng_t *rng, long int nin, double pp) {
  return Rf_rbinom(nin, pp);
}

double splicing_rgamma(splicing_rng_t *rng, double a) {
  return Rf_rgamma(a);
}

#else

/* This is from http://www.getreuer.info/home/randmt, 
   licensed under a BSD license 

   Note that this code was not very well tested, but since we only it
   for generating a random starting point for our Gibbs sampler, 
   it is probably OK. */

double splicing_rgamma(splicing_rng_t *rng, double a) {

  const double b = 1.0;
  const double d = a - 1.0/3.0;
  const double c = 1.0/sqrt(9*d);
  double x, v, u;
    
  if(a >= 1) {

    do {
      do {
	x = splicing_rng_get_normal(rng, /*mean=*/ 0, /*sd=*/ 1);
	v = 1 + c*x;
      } while (v <= 0);
      v = v*v*v;
      u = splicing_rng_get_unif(rng, 0, 1);
      x = x*x;
    } while (u >= 1 - 0.0331*x*x  && log(u) >= x/2 + d*(1 - v + log(v)));

    return b*d*v;

  } else if(a > 0) {

    return splicing_rgamma(rng, 1 + a) * 
      pow(splicing_rng_get_unif(rng, 0, 1), 1/a);

  } else {

    return 0;

  }
}

#endif

/**********************************************************
 * Testing purposes                                       *
 *********************************************************/

/* int main() { */

/*   int i; */

/*   for (i=0; i<1000; i++) { */
/*     printf("%li ", RNG_INTEGER(1,10)); */
/*   } */
/*   printf("\n"); */

/*   for (i=0; i<1000; i++) { */
/*     printf("%f ", RNG_UNIF(0,1)); */
/*   } */
/*   printf("\n"); */

/*   for (i=0; i<1000; i++) { */
/*     printf("%f ", RNG_NORMAL(0,5)); */
/*   } */
/*   printf("\n"); */

/*   return 0; */
/* } */
