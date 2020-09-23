/*
 * malm11.c 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "malm.h"

void
malm11(
  const int *dim,
  const uv *y,
  const double *w,
  uv *E,
  double *alpha,
  double *beta0,
  double *phi0,
  double *L_,
  const int *iopt_,
  const double *dopt_
  )
{
  const int n = dim[0], m = dim[1];
  const int family = iopt_[0];
  const int verbose = iopt_[1];
  const int iter_max = iopt_[2];
  const double epsilon = dopt_[0];

  // init
  for(int i = 0; i < n*m; i++ )
    E[i].u = E[i].v = 1.0;
  for(int i = 0; i < n; i++ )
    alpha[i] = (w[i]==0?NAN:1);
  for(int j = 0; j < m; j++ )
    beta0[j] = phi0[j] = 1;

  int nnz = 0;
  for(int i = 0; i < n; i++ )
    if( w[i] > 0 ) nnz++;
  const double d_adj = 1/(1-(double)(nnz+m-1)/(nnz*m));

  double *Sy = (double*)malloc(sizeof(double)*2*n);
  double *Su = Sy+n;

  for(int ii = 0; ii < 2*n; ii++ )
    Sy[ii] = 0;

  double Lold = -INFINITY;
  for(int iter = 1; iter <= iter_max; iter++ )
    {
    for(int j = 0; j < m; j++ )
      {
      const uv *yj = y + j*n;
      uv *Ej = E + j*n;
      double Syj = 0, Suj = 0;
      for(int i = 0; i < n; i++ )
        {
        if(w[i] == 0) continue;
        double W = w[i] * alpha[i] / Ej[i].v;
        Syj += W * yj[i].u;
        Suj += W * Ej[i].u;
        }
      if( Suj > 0 )
        beta0[j] *= Syj/Suj;
      else
        beta0[j] = 0;
      }
    
    for(int i = 0; i < n; i++ )
      Sy[i] = Su[i] = 0;
    
    for(int j = 0; j < m; j++ )
      {
      const uv *yj = y + j*n;
      uv *Ej = E + j*n;
      for(int i = 0; i < n; i++ )
        {
        if(w[i]==0) continue;
        double W = beta0[j] / Ej[i].v;  // w[i] cancels; omitted
        Sy[i] += W * yj[i].u;
        Su[i] += W * Ej[i].u;
        }
      }

    // constrain sum log alpha = 0
    //
    double Slog = 0; int Sw = 0;
    for(int i = 0; i < n; i++ )
      {
      if( w[i] == 0 ) continue;
      if( Su[i] > 0 )
        alpha[i] *= Sy[i]/Su[i];
      else
        alpha[i] = 0;
      Slog += log(alpha[i]);
      Sw++;
      }
    double geomean_alpha = exp( -Slog/Sw );
    for(int i = 0; i < n; i++)
      {
      if(w[i]==0) continue;
      alpha[i] *= geomean_alpha;
      }

    for(int j = 0; j < m; j++ )
      {
      uv *Ej = E + j*n;
      for(int i = 0; i < n; i++ )
        Ej[i].u = alpha[i]*beta0[j];
      }

    // dispersion
    double L = 0;
    for(int j = 0; j < m; j++ )
      {
      const uv *yj = y + j*n;
      uv *Ej = E + j*n;
      double Sd = 0, Ss = 0;
      for(int i = 0; i < n; i++ )
        {
        if(w[i]==0) continue;
        double d = yj[i].u - Ej[i].u;
        d *= d * d_adj;
        double W = Ej[i].u/Ej[i].v;
        W *= w[i] * W;
        Sd += W * d;
        Ss += W * Ej[i].v;
        }
      if(Ss > 0 )
        phi0[j] *= Sd/Ss;
      else
        phi0[j] = 0;

      if( family == 1 ) // negative binomial
        for(int i = 0; i < n; i++ )
          {
          if( w[i]==0 ) continue;
          Ej[i].v = Ej[i].u + phi0[j] * Ej[i].u * Ej[i].u;
          }
      else // RE gamma
        for(int i = 0; i < n; i++ )
          {
          if( w[i]==0 ) continue;
          Ej[i].v = yj[i].v + phi0[j] * Ej[i].u * Ej[i].u;
          }

      L += PLL( n, w, yj, Ej );
      }

    double change = fabs(L-Lold)/(m*nnz-(m+nnz-1));
    if( verbose == 2 ) mess("[%3d] %g %g\n",iter, L, change);
    if(change < epsilon ) break;
    Lold = L;
    } // iter 

  for(int i = 0; i < n; i++ )
    alpha[i] = log(alpha[i]);
  
  for(int j = 0; j < m; j++ )
    beta0[j] = log(beta0[j] );
  
  if(L_) L_[0] = Lold;
  free(Sy);
}
