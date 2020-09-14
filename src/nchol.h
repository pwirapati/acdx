/*
nchol.h - cholesky decomposition

*/

#ifndef _NCHOL_H_
#define _NCHOL_H_

#include <math.h>


#define NCHOL_AUTOTOL (-1)
#define DOUBLE_NEG_EPS (1.110223e-16) // from R on Intel

static inline double
nchol_autotol(
  const int n,
  const double *x
  )
{
  double maxdiag = -INFINITY;
  const double *xj = x;
  for(int j = 0; j < n; j++, xj += j )
    if( xj[j] > maxdiag ) maxdiag = xj[j];
  return DOUBLE_NEG_EPS * n * maxdiag;
}

static inline void
nchol_add_off(
  const int k,
  const double *R,
  double *x
  )
{
  const double *Ri = R;
  for(int i = 0; i < k-1; i++, Ri += i )
    {
    if( Ri[i] == 0 ) 
      { x[i] = 0; continue; }

    for(int j = 0; j < i; j++ )
      x[i] -= x[j]*Ri[j]; 
    
    x[i] /= Ri[i];
    }
}

static inline void
nchol_add_diag(
  const int k,
  double *x, 
  const double tol
  )
{
  for(int i = 0; i < k-1; i++ )
    x[k-1] -= x[i]*x[i];
  if(x[k-1] > tol )
    x[k-1] = sqrt(x[k-1]);
  else
    x[k-1] = 0;
}

static inline void
nchol_add(
  const int k,
  const double *R,
  double *x,
  const double tol
  )
{
  nchol_add_off(k,R,x);
  nchol_add_diag(k,x,tol);
}

static inline int
nchol(
  const int k,
  double *R,
  double tol
  )
{
  if(tol < 0 ) tol = nchol_autotol(k,R);
  int nskip = 0;
  double *Ri = R;
  for(int i = 0; i < k; i++, Ri += i )
    {
    nchol_add(i+1,R,Ri,tol);
    if(Ri[i]==0) nskip++;
    }
  return nskip;
}

static inline void
nchol_fsub(
  int n,
  const double *R,
  double *y
  )
{
  const double *Rj = R;
  for(int j = 0; j < n; j++, Rj += j )
    if( Rj[j] > 0 )
      {
      for(int k = 0; k < j; k++ )
        y[j] -= y[k]*Rj[k];
      y[j] /= Rj[j];
      }
    else
      y[j] = 0;
}

static inline void
nchol_bsub(
  int n,
  const double *R,
  double *y
  )
{
  const double *Rj = R + n*(n-1)/2;
  for(int j = n-1; j >= 0; Rj -= j, j-- )
    if( Rj[j] > 0 )
      {
      y[j] /= Rj[j];
      for(int k = 0; k < j; k++ )
        y[k] -= y[j]*Rj[k];
      }
    else
      y[j] = 0;
}

static inline double
nchol_ssq(
 const int n,
 const double *x
 )
{
  double ssq = 0;
  for(int i = 0; i < n; i++ )
    ssq += x[i]*x[i];
  return ssq;
}

#endif
