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

typedef struct {
  int i, j;
  long int x[31];
} splicing_i_rng_glibc2_state_t;

unsigned long int splicing_i_rng_glibc2_get(int *i, int *j, int n, 
					  long int *x) {
  long int k;

  x[*i] += x[*j];
  k = (x[*i] >> 1) & 0x7FFFFFFF;
  
  (*i)++;
  if (*i == n) {
    *i = 0;
  }
  
  (*j)++ ;
  if (*j == n) {
    *j = 0;
  }

  return k;
}

unsigned long int splicing_rng_glibc2_get(void *vstate) {
  splicing_i_rng_glibc2_state_t *state = 
    (splicing_i_rng_glibc2_state_t*) vstate;
  return splicing_i_rng_glibc2_get(&state->i, &state->j, 31, state->x);
}

double splicing_rng_glibc2_get_real(void *state) {
  return splicing_rng_glibc2_get(state) / 2147483648.0;
}

/* this function is independent of the bit size */

void splicing_i_rng_glibc2_init(long int *x, int n, 
			      unsigned long int s) {
  int i;
  
  if (s==0) { s=1; }
  
  x[0] = s;
  for (i=1 ; i<n ; i++) {
    const long int h = s / 127773;
    const long int t = 16807 * (s - h * 127773) - h * 2836;
    if (t < 0) {
      s = t + 2147483647 ;
    } else { 
      s = t ;
    }
    
    x[i] = s ;
  }
}

int splicing_rng_glibc2_seed(void *vstate, unsigned long int seed) {
  splicing_i_rng_glibc2_state_t *state = 
    (splicing_i_rng_glibc2_state_t*) vstate;
  int i;
  
  splicing_i_rng_glibc2_init(state->x, 31, seed);
  
  state->i=3;
  state->j=0;
  
  for (i=0;i<10*31; i++) {
    splicing_rng_glibc2_get(state);
  }
  
  return 0;
}

int splicing_rng_glibc2_init(void **state) {
  splicing_i_rng_glibc2_state_t *st;

  st=splicing_Calloc(1, splicing_i_rng_glibc2_state_t);
  if (!st) {
    SPLICING_ERROR("Cannot initialize RNG", SPLICING_ENOMEM);
  }
  (*state)=st;

  splicing_rng_glibc2_seed(st, 0);
  
  return 0;
}

void splicing_rng_glibc2_destroy(void *vstate) {
  splicing_i_rng_glibc2_state_t *state = 
    (splicing_i_rng_glibc2_state_t*) vstate;
  splicing_Free(state);
}

/**
 * \var splicing_rngtype_glibc2
 * \brief The random number generator type introduced in GNU libc 2
 * 
 * It is a linear feedback shift register generator with a 128-byte
 * buffer. This generator was the default prior to splicing version 0.6,
 * at least on systems relying on GNU libc.
 * 
 * This generator was ported from the GNU Scientific Library.
 */

splicing_rng_type_t splicing_rngtype_glibc2 = {
  /* name= */      "LIBC",
  /* min=  */      0,
  /* max=  */      RAND_MAX,
  /* init= */      splicing_rng_glibc2_init,
  /* destroy= */   splicing_rng_glibc2_destroy,
  /* seed= */      splicing_rng_glibc2_seed,
  /* get= */       splicing_rng_glibc2_get,
  /* get_real= */  splicing_rng_glibc2_get_real,
  /* get_norm= */  0,
  /* get_geom= */  0,
  /* get_binom= */ 0,
  /* get_gamma= */ 0
};

/* ------------------------------------ */

typedef struct {
  unsigned long int x;
} splicing_i_rng_rand_state_t;

unsigned long int splicing_rng_rand_get(void *vstate) {
  splicing_i_rng_rand_state_t *state = vstate;
  state->x = (1103515245 * state->x + 12345) & 0x7fffffffUL;
  return state->x;
}

double splicing_rng_rand_get_real(void *vstate) {
  return splicing_rng_rand_get (vstate) / 2147483648.0 ;
}

int splicing_rng_rand_seed(void *vstate, unsigned long int seed) {
  splicing_i_rng_rand_state_t *state = vstate;
  state->x = seed;
  return 0;
}

int splicing_rng_rand_init(void **state) {
  splicing_i_rng_rand_state_t *st;

  st=splicing_Calloc(1, splicing_i_rng_rand_state_t);
  if (!st) {
    SPLICING_ERROR("Cannot initialize RNG", SPLICING_ENOMEM);
  }
  (*state)=st;
  
  splicing_rng_rand_seed(st, 0);
  
  return 0;
}

void splicing_rng_rand_destroy(void *vstate) {
  splicing_i_rng_rand_state_t *state = 
    (splicing_i_rng_rand_state_t*) vstate;
  splicing_Free(state);  
}  

/**
 * \var splicing_rngtype_rand
 * \brief The old BSD rand/stand random number generator
 * 
 * The sequence is 
 *     x_{n+1} = (a x_n + c) mod m 
 * with a = 1103515245, c = 12345 and m = 2^31 = 2147483648. The seed
 * specifies the initial value, x_1.
 * 
 * The theoretical value of x_{10001} is 1910041713.
 *
 *  The period of this generator is 2^31.
 * 
 * This generator is not very good -- the low bits of successive
 * numbers are correlated.
 * 
 * This generator was ported from the GNU Scientific Library.
 */

splicing_rng_type_t splicing_rngtype_rand = {
  /* name= */      "RAND",
  /* min=  */      0,
  /* max=  */      0x7fffffffUL,
  /* init= */      splicing_rng_rand_init,
  /* destroy= */   splicing_rng_rand_destroy,
  /* seed= */      splicing_rng_rand_seed,
  /* get= */       splicing_rng_rand_get,
  /* get_real= */  splicing_rng_rand_get_real,
  /* get_norm= */  0,
  /* get_geom= */  0,
  /* get_binom= */ 0,
  /* get_gamma= */ 0
};

/* ------------------------------------ */

#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;   

/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;   

typedef struct {
  unsigned long mt[N];
  int mti;
} splicing_i_rng_mt19937_state_t;

unsigned long int splicing_rng_mt19937_get(void *vstate) {
  splicing_i_rng_mt19937_state_t *state = vstate;

  unsigned long k ;
  unsigned long int *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

  if (state->mti >= N) {
    /* generate N words at one time */
    int kk;
    
    for (kk = 0; kk < N - M; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
    }
    for (; kk < N - 1; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
    }
    
    {
      unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
    }
    
    state->mti = 0;
  }

#undef MAGIC

  /* Tempering */
  
  k = mt[state->mti];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680UL;
  k ^= (k << 15) & 0xefc60000UL;
  k ^= (k >> 18);

  state->mti++;

  return k;
}

double splicing_rng_mt19937_get_real(void *vstate) {
  return splicing_rng_mt19937_get (vstate) / 4294967296.0 ;
}

int splicing_rng_mt19937_seed(void *vstate, unsigned long int seed) {
  splicing_i_rng_mt19937_state_t *state = vstate;
  int i;

  memset(state, 0, sizeof(splicing_i_rng_mt19937_state_t));

  if (seed == 0) {
    seed = 4357;   /* the default seed is 4357 */
  }
  state->mt[0]= seed & 0xffffffffUL;

  for (i = 1; i < N; i++) {
    /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
       Ed. p.106 for multiplier. */
    state->mt[i] =
      (1812433253UL * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i);
    state->mt[i] &= 0xffffffffUL;
  }
  
  state->mti = i;
  return 0;
}

int splicing_rng_mt19937_init(void **state) {
  splicing_i_rng_mt19937_state_t *st;

  st=splicing_Calloc(1, splicing_i_rng_mt19937_state_t);
  if (!st) {
    SPLICING_ERROR("Cannot initialize RNG", SPLICING_ENOMEM);
  }
  (*state)=st;
  
  splicing_rng_mt19937_seed(st, 0);
  
  return 0;
}

void splicing_rng_mt19937_destroy(void *vstate) {
  splicing_i_rng_mt19937_state_t *state = 
    (splicing_i_rng_mt19937_state_t*) vstate;
  splicing_Free(state);  
}  

/** 
 * \var splicing_rngtype_mt19937
 * \brief The MT19937 random number generator
 * 
 * The MT19937 generator of Makoto Matsumoto and Takuji Nishimura is a
 * variant of the twisted generalized feedback shift-register
 * algorithm, and is known as the “Mersenne Twister” generator. It has
 * a Mersenne prime period of 2^19937 - 1 (about 10^6000) and is
 * equi-distributed in 623 dimensions. It has passed the diehard
 * statistical tests. It uses 624 words of state per generator and is
 * comparable in speed to the other generators. The original generator
 * used a default seed of 4357 and choosing s equal to zero in
 * gsl_rng_set reproduces this. Later versions switched to 5489 as the
 * default seed, you can choose this explicitly via splicing_rng_seed
 * instead if you require it. 
 * 
 * For more information see,
 * Makoto Matsumoto and Takuji Nishimura, “Mersenne Twister: A
 * 623-dimensionally equidistributed uniform pseudorandom number
 * generator”. ACM Transactions on Modeling and Computer Simulation,
 * Vol. 8, No. 1 (Jan. 1998), Pages 3–30 
 * 
 * The generator splicing_rngtype_mt19937 uses the second revision of the
 * seeding procedure published by the two authors above in 2002. The
 * original seeding procedures could cause spurious artifacts for some
 * seed values.
 * 
 * This generator was ported from the GNU Scientific Library.
 */

splicing_rng_type_t splicing_rngtype_mt19937 = {
  /* name= */      "MT19937",
  /* min=  */      0,
  /* max=  */      0xffffffffUL,
  /* init= */      splicing_rng_mt19937_init,
  /* destroy= */   splicing_rng_mt19937_destroy,
  /* seed= */      splicing_rng_mt19937_seed,
  /* get= */       splicing_rng_mt19937_get,
  /* get_real= */  splicing_rng_mt19937_get_real,
  /* get_norm= */  0,
  /* get_geom= */  0,
  /* get_binom= */ 0
};

#undef N
#undef M

/* ------------------------------------ */

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

#ifndef USING_R

#define addr(a) (&a)

splicing_i_rng_mt19937_state_t splicing_i_rng_default_state;

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
  addr(splicing_rngtype_mt19937),
  addr(splicing_i_rng_default_state),
  /* def= */ 1
};

#undef addr

#endif

/* ------------------------------------ */

#ifdef USING_R

double  unif_rand(void);
double  norm_rand(void);
double  Rf_rgeom(double);
double  Rf_rbinom(double, double);
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

double splicing_rng_R_get_geom(void *state, double p) {
  return Rf_rgeom(p);
}
 
double splicing_rng_R_get_binom(void *state, long int n,
				     double p) {
  return Rf_rbinom(n, p);
}

double splicing_rng_R_get_gamma(void *state, double a, double scale) {
  return Rf_rgamma(a, scale);
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
  /* get_geom= */  splicing_rng_R_get_geom,
  /* get_binom= */ splicing_rng_R_get_binom,
  /* get_gamma= */ splicing_rng_R_get_gamma
};

splicing_rng_t splicing_rng_default = { 
  &splicing_rngtype_R,
  0,
  /* def= */ 1
};

#endif

/* ------------------------------------ */

double splicing_norm_rand(splicing_rng_t *rng);
double splicing_rgeom(splicing_rng_t *rng, double p);
double splicing_rbinom(splicing_rng_t *rng, long int nin, double pp);
double splicing_rgamma(splicing_rng_t *rng, double a, double scale);

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
  } else {
    return splicing_norm_rand(rng)*s+m;
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

/** 
 * \function splicing_rng_get_geom
 * Generate geometrically distributed random numbers
 * 
 * \param rng Pointer to the RNG to use. Use \ref splicing_rng_default
 *        here to use the default splicing RNG.
 * \param p The probability of success in each trial. Must be larger
 *        than zero and smaller or equal to 1.
 * \return The generated geometrically distributed random number.
 * 
 * Time complexity: depends on the type of the RNG.
 */

double splicing_rng_get_geom(splicing_rng_t *rng, double p) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_geom) {
    return type->get_geom(rng->state, p);
  } else {
    return splicing_rgeom(rng, p);
  }
}

/** 
 * \function splicing_rng_get_binom
 * Generate binomially distributed random numbers
 * 
 * \param rng Pointer to the RNG to use. Use \ref splicing_rng_default
 *        here to use the default splicing RNG.
 * \param n Number of observations.
 * \param p Probability of an event.
 * \return The generated binomially distributed random number.
 * 
 * Time complexity: depends on the type of the RNG.
 */

double splicing_rng_get_binom(splicing_rng_t *rng, long int n, 
				   double p) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_binom) {
    return type->get_binom(rng->state, n, p);
  } else {
    return splicing_rbinom(rng, n, p);
  }
}

double splicing_rng_get_gamma(splicing_rng_t *rng, double a, 
			      double scale) {
  const splicing_rng_type_t *type=rng->type;
  if (type->get_gamma) {
    return type->get_gamma(rng->state, a, scale);
  } else {
    return splicing_rgamma(rng, a, scale);
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

double splicing_rgamma(splicing_rng_t *rng, double a, double scale) {
  return Rf_rgamma(a, scale);
}

#else

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
 *  based on AS 111 (C) 1977 Royal Statistical Society
 *  and   on AS 241 (C) 1988 Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *	double qnorm5(double p, double mu, double sigma,
 *		      int lower_tail, int log_p)
 *            {qnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *	Compute the quantile function for the normal distribution.
 *
 *	For small to moderate probabilities, algorithm referenced
 *	below is used to obtain an initial approximation which is
 *	polished with a final Newton step.
 *
 *	For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *	Beasley, J. D. and S. G. Springer (1977).
 *	Algorithm AS 111: The percentage points of the normal distribution,
 *	Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2004  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifdef _MSC_VER
#  define ML_POSINF SPLICING_INFINITY
#  define ML_NEGINF -SPLICING_INFINITY
#  define ML_NAN    SPLICING_NAN
#else
#  define ML_POSINF	(1.0 / 0.0)
#  define ML_NEGINF	((-1.0) / 0.0)
#  define ML_NAN		(0.0 / 0.0)
#endif

#define ML_ERROR(x)	/* nothing */
#define ML_UNDERFLOW	(DBL_MIN * DBL_MIN)
#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occurred (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

/* Wilcoxon Signed Rank Distribution */

#define SIGNRANK_MAX 50

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0       	Rf_bd0
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define i1mach		Rf_i1mach
#define gammalims	Rf_gammalims
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define stirlerr       	Rf_stirlerr

	/* Chebyshev Series */

int	chebyshev_init(double*, int, double);
double	chebyshev_eval(double, const double *, const int);

	/* Gamma and Related Functions */

void	gammalims(double*, double*);
double	lgammacor(double); /* log(gamma) correction */
double  stirlerr(double);  /* Stirling expansion "error" */

double	lfastchoose(double, double);

double  bd0(double, double);

/* Consider adding these two to the API (Rmath.h): */
double	dbinom_raw(double, double, double, double, int);
double	dpois_raw (double, double, int);
double  pnchisq_raw(double, double, double, double, double, int);

int	i1mach(int);

/* From toms708.c */
void bratio(double a, double b, double x, double y, 
	    double *w, double *w1, int *ierr);


#endif /* MATHLIB_PRIVATE_H */


	/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

#define R_D_Lval(p)	(lower_tail ? (p) : (1 - (p)))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (1 - (p)) : (p))	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (1 - (p)))/* [log](1-p) */

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

/*till 1.8.x:
 * #define R_DT_val(x)	R_D_val(R_D_Lval(x))
 * #define R_DT_Cval(x)	R_D_val(R_D_Cval(x)) */
#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)	(lower_tail ? R_D_Clog(x) : R_D_val(x))

/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))

#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/
#define R_DT_Log(p)	(lower_tail? (p) : R_Log1_Exp(p))
/* ==   R_DT_log when we already "know" log_p == TRUE :*/

#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_ERR_return_NAN

/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x) 	  (fabs((x) - floor((x)+0.5)) > 1e-7)
/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))

#define R_D_nonint_check(x) 				\
   if(R_D_nonint(x)) {					\
	MATHLIB_WARNING("non-integer x = %f", x);	\
	return R_D__0;					\
   }

double splicing_qnorm5(double p, double mu, double sigma, int lower_tail, int log_p)
{
    double p_, q, r, val;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
	return p + mu + sigma;
#endif
    if (p == R_DT_0)	return ML_NEGINF;
    if (p == R_DT_1)	return ML_POSINF;
    R_Q_P01_check(p);

    if(sigma  < 0)	ML_ERR_return_NAN;
    if(sigma == 0)	return mu;

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

	/* r = min(p, 1-p) < 0.075 */
	if (q > 0)
	    r = R_DT_CIv(p);/* 1-p */
	else
	    r = p_;/* = R_DT_Iv(p) ^=  p */

	r = sqrt(- ((log_p &&
		     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
		    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

	if(q < 0.0)
	    val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

double fsign(double x, double y)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(y))
        return x + y;
#endif
    return ((y >= 0) ? fabs(x) : -fabs(x));
}

int imax2(int x, int y)
{
    return (x < y) ? y : x;
}

int imin2(int x, int y)
{
    return (x < y) ? x : y;
}

#ifdef HAVE_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#elif HAVE_WORKING_FINITE
/* include header needed to define finite() */
#  ifdef HAVE_IEEE754_H
#   include <ieee754.h>         /* newer Linuxen */
#  else
#   ifdef HAVE_IEEEFP_H
#    include <ieeefp.h>         /* others [Solaris], .. */
#   endif
#  endif
# define R_FINITE(x)    finite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif

int R_finite(double x)
{
#ifdef HAVE_ISFINITE
    return isfinite(x);
#elif HAVE_WORKING_FINITE
    return finite(x);
#else
/* neither finite nor isfinite work. Do we really need the AIX exception? */
# ifdef _AIX
#  include <fp.h>
     return FINITE(x);
# elif defined(_MSC_VER)
     return _finite(x);
#else
    return (!isnan(x) & (x != 1/0.0) & (x != -1.0/0.0));
# endif
#endif
}

int R_isnancpp(double x)
{
   return (isnan(x)!=0);
}

#ifdef __cplusplus
  int R_isnancpp(double); /* in arithmetic.c */
#  define ISNAN(x)     R_isnancpp(x)
#else
#  define ISNAN(x)     (isnan(x)!=0)
#endif

double splicing_norm_rand(splicing_rng_t *rng) {
  
  double u1;

#define BIG 134217728 /* 2^27 */
  /* unif_rand() alone is not of high enough precision */
  u1 = splicing_rng_get_unif01(rng);
  u1 = (int)(BIG*u1) + splicing_rng_get_unif01(rng);
  return splicing_qnorm5(u1/BIG, 0.0, 1.0, 1, 0);
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2002 the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double exp_rand(void);
 *
 *  DESCRIPTION
 *
 *    Random variates from the standard exponential distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1972).
 *    Computer methods for sampling from the exponential and
 *    normal distributions.
 *    Comm. ACM, 15, 873-882.
 */

double splicing_exp_rand(splicing_rng_t *rng)
{
    /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
    /* The highest n (here 8) is determined by q[n-1] = 1.0 */
    /* within standard precision */
    const double q[] =
    {
	0.6931471805599453,
	0.9333736875190459,
	0.9888777961838675,
	0.9984959252914960,
	0.9998292811061389,
	0.9999833164100727,
	0.9999985691438767,
	0.9999998906925558,
	0.9999999924734159,
	0.9999999995283275,
	0.9999999999728814,
	0.9999999999985598,
	0.9999999999999289,
	0.9999999999999968,
	0.9999999999999999,
	1.0000000000000000
    };
    double a, u, ustar, umin;
    int i;

    a = 0.;
    /* precaution if u = 0 is ever returned */
    u = splicing_rng_get_unif01(rng);
    while(u <= 0.0 || u >= 1.0) u = splicing_rng_get_unif01(rng);
    for (;;) {
	u += u;
	if (u > 1.0)
	    break;
	a += q[0];
    }
    u -= 1.;

    if (u <= q[0])
	return a + u;

    i = 0;
    ustar = splicing_rng_get_unif01(rng);
    umin = ustar;
    do {
        ustar = splicing_rng_get_unif01(rng);
	if (ustar < umin)
	    umin = ustar;
	i++;
    } while (u > q[i]);
    return a + umin * q[0];
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2001 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rpois(double lambda)
 *
 *  DESCRIPTION
 *
 *    Random variates from the Poisson distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1982).
 *    Computer generation of Poisson deviates
 *    from modified normal distributions.
 *    ACM Trans. Math. Software 8, 163-179.
 */

#define a0	-0.5
#define a1	 0.3333333
#define a2	-0.2500068
#define a3	 0.2000118
#define a4	-0.1661269
#define a5	 0.1421878
#define a6	-0.1384794
#define a7	 0.1250060

#define one_7	0.1428571428571428571
#define one_12	0.0833333333333333333
#define one_24	0.0416666666666666667

#define repeat for(;;)

#define FALSE 0
#define TRUE  1
#define M_1_SQRT_2PI    0.398942280401432677939946059934     /* 1/sqrt(2pi) */

double splicing_rpois(splicing_rng_t *rng, double mu)
{
    /* Factorial Table (0:9)! */
    const double fact[10] =
    {
	1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };

    /* These are static --- persistent between calls for same mu : */
    static int l, m;

    static double b1, b2, c, c0, c1, c2, c3;
    static double pp[36], p0, p, q, s, d, omega;
    static double big_l;/* integer "w/o overflow" */
    static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/

    /* Local Vars  [initialize some for -Wall]: */
    double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
    double pois = -1.;
    int k, kflag, big_mu, new_big_mu = FALSE;

    if (!R_FINITE(mu))
	ML_ERR_return_NAN;

    if (mu <= 0.)
	return 0.;

    big_mu = mu >= 10.;
    if(big_mu)
	new_big_mu = FALSE;

    if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

	if (big_mu) {
	    new_big_mu = TRUE;
	    /* Case A. (recalculation of s,d,l	because mu has changed):
	     * The Poisson probabilities pk exceed the discrete normal
	     * probabilities fk whenever k >= m(mu).
	     */
	    muprev = mu;
	    s = sqrt(mu);
	    d = 6. * mu * mu;
	    big_l = floor(mu - 1.1484);
	    /* = an upper bound to m(mu) for all mu >= 10.*/
	}
	else { /* Small mu ( < 10) -- not using normal approx. */

	    /* Case B. (start new table and calculate p0 if necessary) */

	    /*muprev = 0.;-* such that next time, mu != muprev ..*/
	    if (mu != muprev) {
		muprev = mu;
		m = imax2(1, (int) mu);
		l = 0; /* pp[] is already ok up to pp[l] */
		q = p0 = p = exp(-mu);
	    }

	    repeat {
		/* Step U. uniform sample for inversion method */
	        u = splicing_rng_get_unif01(rng);
		if (u <= p0)
		    return 0.;

		/* Step T. table comparison until the end pp[l] of the
		   pp-table of cumulative Poisson probabilities
		   (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
		if (l != 0) {
		    for (k = (u <= 0.458) ? 1 : imin2(l, m);  k <= l; k++)
			if (u <= pp[k])
			    return (double)k;
		    if (l == 35) /* u > pp[35] */
			continue;
		}
		/* Step C. creation of new Poisson
		   probabilities p[l..] and their cumulatives q =: pp[k] */
		l++;
		for (k = l; k <= 35; k++) {
		    p *= mu / k;
		    q += p;
		    pp[k] = q;
		    if (u <= q) {
			l = k;
			return (double)k;
		    }
		}
		l = 35;
	    } /* end(repeat) */
	}/* mu < 10 */

    } /* end {initialize persistent vars} */

/* Only if mu >= 10 : ----------------------- */

    /* Step N. normal sample */
    g = mu + s * splicing_norm_rand(rng);/* norm_rand() ~ N(0,1), standard normal */

    if (g >= 0.) {
	pois = floor(g);
	/* Step I. immediate acceptance if pois is large enough */
	if (pois >= big_l)
	    return pois;
	/* Step S. squeeze acceptance */
	fk = pois;
	difmuk = mu - fk;
	u = splicing_rng_get_unif01(rng); /* ~ U(0,1) - sample */
	if (d * u >= difmuk * difmuk * difmuk)
	    return pois;
    }

    /* Step P. preparations for steps Q and H.
       (recalculations of parameters if necessary) */

    if (new_big_mu || mu != muprev2) {
        /* Careful! muprev2 is not always == muprev
	   because one might have exited in step I or S
	   */
        muprev2 = mu;
	omega = M_1_SQRT_2PI / s;
	/* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
	 * approximations to the discrete normal probabilities fk. */

	b1 = one_24 / mu;
	b2 = 0.3 * b1 * b1;
	c3 = one_7 * b1 * b2;
	c2 = b2 - 15. * c3;
	c1 = b1 - 6. * b2 + 45. * c3;
	c0 = 1. - b1 + 3. * b2 - 15. * c3;
	c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
    }

    if (g >= 0.) {
	/* 'Subroutine' F is called (kflag=0 for correct return) */
	kflag = 0;
	goto Step_F;
    }


    repeat {
	/* Step E. Exponential Sample */

	E = splicing_exp_rand(rng);/* ~ Exp(1) (standard exponential) */

	/*  sample t from the laplace 'hat'
	    (if t <= -0.6744 then pk < fk for all mu >= 10.) */
	u = 2 * splicing_rng_get_unif01(rng) - 1.;
	t = 1.8 + fsign(E, u);
	if (t > -0.6744) {
	    pois = floor(mu + s * t);
	    fk = pois;
	    difmuk = mu - fk;

	    /* 'subroutine' F is called (kflag=1 for correct return) */
	    kflag = 1;

	  Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

	    if (pois < 10) { /* use factorials from table fact[] */
		px = -mu;
		py = pow(mu, pois) / fact[(int)pois];
	    }
	    else {
		/* Case pois >= 10 uses polynomial approximation
		   a0-a7 for accuracy when advisable */
		del = one_12 / fk;
		del = del * (1. - 4.8 * del * del);
		v = difmuk / fk;
		if (fabs(v) <= 0.25)
		    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
					  v + a3) * v + a2) * v + a1) * v + a0)
			- del;
		else /* |v| > 1/4 */
		    px = fk * log(1. + v) - difmuk - del;
		py = M_1_SQRT_2PI / sqrt(fk);
	    }
	    x = (0.5 - difmuk) / s;
	    x *= x;/* x^2 */
	    fx = -0.5 * x;
	    fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
	    if (kflag > 0) {
		/* Step H. Hat acceptance (E is repeated on rejection) */
		if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
		    break;
	    } else
		/* Step Q. Quotient acceptance (rare case) */
		if (fy - u * fy <= py * exp(px - fx))
		    break;
	}/* t > -.67.. */
    }
    return pois;
}

#undef a0
#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef a7

double splicing_rgeom(splicing_rng_t *rng, double p) {
    if (ISNAN(p) || p <= 0 || p > 1) ML_ERR_return_NAN;

    return splicing_rpois(rng, splicing_exp_rand(rng) * ((1 - p) / p));
}

/* This is from nmath/rbinom.c */

double splicing_rbinom(splicing_rng_t *rng, long int nin, double pp)
{
    /* FIXME: These should become THREAD_specific globals : */

    static double c, fm, npq, p1, p2, p3, p4, qn;
    static double xl, xll, xlr, xm, xr;

    static double psave = -1.0;
    static int nsave = -1;
    static int m;

    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    int i,ix,k, n;

    if (!R_FINITE(nin)) ML_ERR_return_NAN;
    n = floor(nin + 0.5);
    if (n != nin) ML_ERR_return_NAN;

    if (!R_FINITE(pp) ||
	/* n=0, p=0, p=1 are not errors <TSL>*/
	n < 0 || pp < 0. || pp > 1.)	ML_ERR_return_NAN;

    if (n == 0 || pp == 0.) return 0;
    if (pp == 1.) return n;

    p = fmin(pp, 1. - pp);
    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);

    /* Setup, perform only when parameters change [using static (globals): */

    /* FIXING: Want this thread safe
       -- use as little (thread globals) as possible
    */
    if (pp != psave || n != nsave) {
	psave = pp;
	nsave = n;
	if (np < 30.0) {
	    /* inverse cdf logic for mean less than 30 */
	    qn = pow(q, (double) n);
	    goto L_np_small;
	} else {
	    ffm = np + p;
	    m = ffm;
	    fm = m;
	    npq = np * q;
	    p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
	    xm = fm + 0.5;
	    xl = xm - p1;
	    xr = xm + p1;
	    c = 0.134 + 20.5 / (15.3 + fm);
	    al = (ffm - xl) / (ffm - xl * p);
	    xll = al * (1.0 + 0.5 * al);
	    al = (xr - ffm) / (xr * q);
	    xlr = al * (1.0 + 0.5 * al);
	    p2 = p1 * (1.0 + c + c);
	    p3 = p2 + c / xll;
	    p4 = p3 + c / xlr;
	}
    } else if (n == nsave) {
	if (np < 30.0)
	    goto L_np_small;
    }

    /*-------------------------- np = n*p >= 30 : ------------------- */
    repeat {
      u = splicing_rng_get_unif01(rng) * p4;
      v = splicing_rng_get_unif01(rng);
      /* triangular region */
      if (u <= p1) {
	  ix = xm - p1 * v + u;
	  goto finis;
      }
      /* parallelogram region */
      if (u <= p2) {
	  x = xl + (u - p1) / c;
	  v = v * c + 1.0 - fabs(xm - x) / p1;
	  if (v > 1.0 || v <= 0.)
	      continue;
	  ix = x;
      } else {
	  if (u > p3) {	/* right tail */
	      ix = xr - log(v) / xlr;
	      if (ix > n)
		  continue;
	      v = v * (u - p3) * xlr;
	  } else {/* left tail */
	      ix = xl + log(v) / xll;
	      if (ix < 0)
		  continue;
	      v = v * (u - p2) * xll;
	  }
      }
      /* determine appropriate way to perform accept/reject test */
      k = abs(ix - m);
      if (k <= 20 || k >= npq / 2 - 1) {
	  /* explicit evaluation */
	  f = 1.0;
	  if (m < ix) {
	      for (i = m + 1; i <= ix; i++)
		  f *= (g / i - r);
	  } else if (m != ix) {
	      for (i = ix + 1; i <= m; i++)
		  f /= (g / i - r);
	  }
	  if (v <= f)
	      goto finis;
      } else {
	  /* squeezing using upper and lower bounds on log(f(x)) */
	  amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
	  ynorm = -k * k / (2.0 * npq);
	  alv = log(v);
	  if (alv < ynorm - amaxp)
	      goto finis;
	  if (alv <= ynorm + amaxp) {
	      /* Stirling's formula to machine accuracy */
	      /* for the final acceptance/rejection test */
	      x1 = ix + 1;
	      f1 = fm + 1.0;
	      z = n + 1 - fm;
	      w = n - ix + 1.0;
	      z2 = z * z;
	      x2 = x1 * x1;
	      f2 = f1 * f1;
	      w2 = w * w;
	      if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
		  goto finis;
	  }
      }
  }

 L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */

  repeat {
     ix = 0;
     f = qn;
     u = splicing_rng_get_unif01(rng);
     repeat {
	 if (u < f)
	     goto finis;
	 if (ix > 110)
	     break;
	 u -= f;
	 ix++;
	 f *= (g / ix - r);
     }
  }
 finis:
    if (psave > 0.5)
	 ix = n - ix;
  return (double)ix;
}

/* From nmath/rgamma.c */

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000--2008 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rgamma(double a, double scale);
 *
 *  DESCRIPTION
 *
 *    Random variates from the gamma distribution.
 *
 *  REFERENCES
 *
 *    [1] Shape parameter a >= 1.  Algorithm GD in:
 *
 *	  Ahrens, J.H. and Dieter, U. (1982).
 *	  Generating gamma variates by a modified
 *	  rejection technique.
 *	  Comm. ACM, 25, 47-54.
 *
 *
 *    [2] Shape parameter 0 < a < 1. Algorithm GS in:
 *
 *	  Ahrens, J.H. and Dieter, U. (1974).
 *	  Computer methods for sampling from gamma, beta,
 *	  poisson and binomial distributions.
 *	  Computing, 12, 223-246.
 *
 *    Input: a = parameter (mean) of the standard gamma distribution.
 *    Output: a variate from the gamma(a)-distribution
 */

double splicing_rgamma(splicing_rng_t *rng, double a, double scale)
{
/* Constants : */
    const static double sqrt32 = 5.656854;
    const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

    /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
     * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
     * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
     */
    const static double q1 = 0.04166669;
    const static double q2 = 0.02083148;
    const static double q3 = 0.00801191;
    const static double q4 = 0.00144121;
    const static double q5 = -7.388e-5;
    const static double q6 = 2.4511e-4;
    const static double q7 = 2.424e-4;

    const static double a1 = 0.3333333;
    const static double a2 = -0.250003;
    const static double a3 = 0.2000062;
    const static double a4 = -0.1662921;
    const static double a5 = 0.1423657;
    const static double a6 = -0.1367177;
    const static double a7 = 0.1233795;

    /* State variables [FIXME for threading!] :*/
    static double aa = 0.;
    static double aaa = 0.;
    static double s, s2, d;    /* no. 1 (step 1) */
    static double q0, b, si, c;/* no. 2 (step 4) */

    double e, p, q, r, t, u, v, w, x, ret_val;

    if (!R_FINITE(a) || !R_FINITE(scale) || a < 0.0 || scale <= 0.0) {
	if(scale == 0.) return 0.;
	ML_ERR_return_NAN;
    }

    if (a < 1.) { /* GS algorithm for parameters a < 1 */
	if(a == 0)
	    return 0.;
	e = 1.0 + exp_m1 * a;
	repeat {
	    p = e * splicing_rng_get_unif01(rng);
	    if (p >= 1.0) {
		x = -log((e - p) / a);
		if (splicing_exp_rand(rng) >= (1.0 - a) * log(x))
		    break;
	    } else {
		x = exp(log(p) / a);
		if (splicing_exp_rand(rng) >= x)
		    break;
	    }
	}
	return scale * x;
    }

    /* --- a >= 1 : GD algorithm --- */

    /* Step 1: Recalculations of s2, s, d if a has changed */
    if (a != aa) {
	aa = a;
	s2 = a - 0.5;
	s = sqrt(s2);
	d = sqrt32 - s * 12.0;
    }
    /* Step 2: t = standard normal deviate,
               x = (s,1/2) -normal deviate. */

    /* immediate acceptance (i) */
    t = splicing_norm_rand(rng);
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0)
	return scale * ret_val;

    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = splicing_rng_get_unif01(rng);
    if (d * u <= t * t * t)
	return scale * ret_val;

    /* Step 4: recalculations of q0, b, si, c if necessary */

    if (a != aaa) {
	aaa = a;
	r = 1.0 / a;
	q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
	       + q2) * r + q1) * r;

	/* Approximation depending on size of parameter a */
	/* The constants in the expressions for b, si and c */
	/* were established by numerical experiments */

	if (a <= 3.686) {
	    b = 0.463 + s + 0.178 * s2;
	    si = 1.235;
	    c = 0.195 / s - 0.079 + 0.16 * s;
	} else if (a <= 13.022) {
	    b = 1.654 + 0.0076 * s2;
	    si = 1.68 / s + 0.275;
	    c = 0.062 / s + 0.024;
	} else {
	    b = 1.77;
	    si = 0.75;
	    c = 0.1515 / s;
	}
    }
    /* Step 5: no quotient test if x not positive */

    if (x > 0.0) {
	/* Step 6: calculation of v and quotient q */
	v = t / (s + s);
	if (fabs(v) <= 0.25)
	    q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
				      + a3) * v + a2) * v + a1) * v;
	else
	    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);


	/* Step 7: quotient acceptance (q) */
	if (log(1.0 - u) <= q)
	    return scale * ret_val;
    }

    repeat {
	/* Step 8: e = standard exponential deviate
	 *	u =  0,1 -uniform deviate
	 *	t = (b,si)-double exponential (laplace) sample */
	e = splicing_exp_rand(rng);
	u = splicing_rng_get_unif01(rng);
	u = u + u - 1.0;
	if (u < 0.0)
	    t = b - si * e;
	else
	    t = b + si * e;
	/* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
	if (t >= -0.71874483771719) {
	    /* Step 10:	 calculation of v and quotient q */
	    v = t / (s + s);
	    if (fabs(v) <= 0.25)
		q = q0 + 0.5 * t * t *
		    ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
		      + a2) * v + a1) * v;
	    else
		q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
	    /* Step 11:	 hat acceptance (h) */
	    /* (if q not positive go to step 8) */
	    if (q > 0.0) {
		w = expm1(q);
		/*  ^^^^^ original code had approximation with rel.err < 2e-7 */
		/* if t is rejected sample again at step 8 */
		if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
		    break;
	    }
	}
    } /* repeat .. until  `t' is accepted */
    x = s + 0.5 * t;
    return scale * x * x;
}

double splicing_rng_get_invchi2(splicing_rng_t *rng, double nu, double tau2) {

  /* This is the same as an Inv-Gamma(nu/2, nu*tau2/2) */
  double alpha = nu / 2.0;
  double beta  = nu * tau2 / 2.0;

  /* And that is the same as 1 / Gamma(alpha, 1/beta) */
  return 1.0 / RNG_GAMMA(alpha, 1.0 / beta);
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
