/* -*- mode: C -*-  */

#ifndef SPLICING_RANDOM_H
#define SPLICING_RANDOM_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#include <stdlib.h>
#include <time.h>

/* The new RNG interface is (somewhat) modelled based on the GSL */

typedef struct splicing_rng_type_t {
  const char *name;
  unsigned long int min;
  unsigned long int max;
  int (*init)(void **state);
  void (*destroy)(void *state);
  int (*seed)(void *state, unsigned long int seed);  
  unsigned long int (*get)(void *state);
  double (*get_real)(void *state);
  double (*get_norm)(void *state);
  double (*get_geom)(void *state, double p);
  double (*get_binom)(void *state, long int n, double p);
  double (*get_gamma)(void *state, double a, double scale);
} splicing_rng_type_t;

typedef struct splicing_rng_t {
  const splicing_rng_type_t *type;
  void *state;
  int def;
} splicing_rng_t;

/* --------------------------------- */

int splicing_rng_init(splicing_rng_t *rng, const splicing_rng_type_t *type);
void splicing_rng_destroy(splicing_rng_t *rng);

int splicing_rng_seed(splicing_rng_t *rng, unsigned long int seed);
unsigned long int splicing_rng_max(splicing_rng_t *rng);
unsigned long int splicing_rng_min(splicing_rng_t *rng);
const char *splicing_rng_name(splicing_rng_t *rng);

long int splicing_rng_get_integer(splicing_rng_t *rng,
				long int l, long int h);
double splicing_rng_get_normal(splicing_rng_t *rng, 
				    double m, double s);
double splicing_rng_get_unif(splicing_rng_t *rng, 
				  double l, double h);
double splicing_rng_get_unif01(splicing_rng_t *rng);
double splicing_rng_get_geom(splicing_rng_t *rng, double p);
double splicing_rng_get_binom(splicing_rng_t *rng, long int n, 
				   double p);
double splicing_rng_get_gamma(splicing_rng_t *rng, double a, 
			      double scale);
unsigned long int splicing_rng_get_int31(splicing_rng_t *rng);

/* --------------------------------- */

extern splicing_rng_type_t splicing_rngtype_glibc2;
extern splicing_rng_type_t splicing_rngtype_rand;
extern splicing_rng_type_t splicing_rngtype_mt19937;
extern splicing_rng_t splicing_rng_default;

void splicing_rng_set_default(splicing_rng_t *rng);

/* --------------------------------- */

#ifdef USING_R

void GetRNGstate(void);
void PutRNGstate(void);
#define RNG_BEGIN()    GetRNGstate()
#define RNG_END()      PutRNGstate()

#else 

#define RNG_BEGIN()      if (splicing_rng_default.def==1) { \
  splicing_rng_seed(&splicing_rng_default, time(0)); \
  splicing_rng_default.def=2; \
  }
#define RNG_END()		/* do nothing */

#endif

#define RNG_INTEGER(l,h) (splicing_rng_get_integer(&splicing_rng_default,(l),(h)))
#define RNG_NORMAL(m,s)  (splicing_rng_get_normal(&splicing_rng_default,(m),(s)))
#define RNG_UNIF(l,h)    (splicing_rng_get_unif(&splicing_rng_default,(l),(h)))
#define RNG_UNIF01()     (splicing_rng_get_unif01(&splicing_rng_default))
#define RNG_GEOM(p)      (splicing_rng_get_geom(&splicing_rng_default,(p)))
#define RNG_BINOM(n,p)   (splicing_rng_get_binom(&splicing_rng_default,(n),(p)))
#define RNG_GAMMA(a,s)   (splicing_rng_get_gamma(&splicing_rng_default,(a),(s)))
#define RNG_INT31()      (splicing_rng_get_int31(&splicing_rng_default))

__END_DECLS

#endif
