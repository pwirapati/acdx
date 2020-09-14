#ifndef _MALM_H_
#define _MALM_H_

#ifndef SHELLPROG 
  #include <R.h> // for Rprintf
  #define mess(...) Rprintf(__VA_ARGS__)
#else
  #define mess(...) fprintf(stderr,__VA_ARGS__)
#endif

typedef struct {
  double u;
  double v;
  }
uv;


static inline void
dset(int n, double *x, double c)
{
  for(int i = 0; i < n; i++ )
    x[i] = c;
}

static inline void
dcp(int n, const double *from, double *to)
{
  for(int i = 0; i < n; i++ )
    to[i] = from[i];
}

static inline void
daccum(int n, double *y, const double *x)
{
  for(int i = 0; i < n; i++ )
    y[i] += x[i];
}

static int
nlevels (int n_ac, const int *x)
{
  int k = -1;
  for(int i = 0; i < n_ac; i++ )
    if( x[i] >= 0 && x[i] > k ) k = x[i];
  return k + 1;
}


// pseudo-log-likelihood
extern double
PLL ( 
  const int n,
  const double *w,
  const uv *Y,
  const uv *E
  )
;

extern void
malm (
  const int *dims_,
  const uv *Y,
  const double *w,
  const double *Xt,
  const int *Gi,
  uv *Em,
  const double *offset,
  double *Beta,
  double *Phi,
  double *L,
  const int *iopt_,
  const double *dopt_,
  int *n_iter
  )
;

extern void
malm11(
  const int *dim,
  const uv *y,
  const double *w,
  uv *E,
  double *alpha,
  double *gamma,
  double *phi,
  double *L_,
  const int *iopt_,
  const double *dopt_
  )
;

#endif  // _MALM_H_
