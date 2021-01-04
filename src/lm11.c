/*
 * lm11.c - normal model with identity link
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "malm.h"

void
lm11(
  const int *dim,
  const uv *y,
  const double *w,
  uv *E,
  double *alpha,
  double *beta0,
  double *tau2,
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
    { E[i].u = 0;  E[i].v = 1; }
  for(int i = 0; i < n; i++ )
    alpha[i] = (w[i]==0?NAN:0);
  for(int j = 0; j < m; j++ )
    { beta0[j] = 0; tau2[j] = 1; }

  int nnz = 0;
  for(int i = 0; i < n; i++ )
    if( w[i] > 0 ) nnz++;
  const double d_adj = 1/(1-(double)(nnz+m-1)/(nnz*m));

  double *Sy = (double*)malloc(sizeof(double)*2*n);
  double *Sw = Sy+n;

  for(int ii = 0; ii < 2*n; ii++ )
    Sy[ii] = 0;

  double Lold = -INFINITY;
  for(int iter = 1; iter <= iter_max; iter++ )
    {
    // update beta0
    for(int j = 0; j < m; j++ )
      {
      const uv *yj = y + j*n;
      uv *Ej = E + j*n;
      double Syj = 0, Swj = 0;
      for(int i = 0; i < n; i++ )
        {
        if(w[i] == 0) continue;
        double W = w[i] / Ej[i].v;
        Syj += W * (yj[i].u - Ej[i].u);
        Swj += W;
        }
      if( Swj > 0 )
        beta0[j] += Syj/Swj;
      else
        beta0[j] = 0;
      }
    
    // update alpha
    for(int i = 0; i < n; i++ )
      Sy[i] = Sw[i] = 0;
    
    for(int j = 0; j < m; j++ )
      {
      const uv *yj = y + j*n;
      uv *Ej = E + j*n;
      for(int i = 0; i < n; i++ )
        {
        if(w[i]==0) continue;
        double W = 1 / Ej[i].v;  // w[i] cancels; omitted
        Sy[i] += W * (yj[i].u - Ej[i].u);
        Sw[i] += W;
        }
      }

    // constrain sum alpha = 0
    //
    double Sa = 0; int Swa = 0;
    for(int i = 0; i < n; i++ )
      {
      if( w[i] == 0 ) continue;
      if( Sw[i] > 0 )
        alpha[i] += Sy[i]/Sw[i];
      else
        alpha[i] += 0;
      Sa += alpha[i];
      Swa++;
      }
    double mean_alpha = Sa/Swa;
    for(int i = 0; i < n; i++)
      {
      if(w[i]==0) continue;
      alpha[i] -= mean_alpha;
      }

    // update mean
    for(int j = 0; j < m; j++ )
      {
      uv *Ej = E + j*n;
      for(int i = 0; i < n; i++ )
        Ej[i].u = alpha[i] + beta0[j];
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
        double W = w[i]/Ej[i].v;
        Sd += W * d /Ej[i].v;
        Ss += W;
        }
      if(Ss > 0 )
        tau2[j] *= Sd/Ss;
      else
        tau2[j] = 0;

      for(int i = 0; i < n; i++ )
        {
        if( w[i]==0 ) continue;
        Ej[i].v = yj[i].v + tau2[j];
        }

      L += PLL( n, w, yj, Ej );
      }

    double change = fabs(L-Lold)/(m*nnz-(m+nnz-1));
    if( verbose == 2 ) mess("[%3d] %g %g\n",iter, L, change);
    if(change < epsilon ) break;
    Lold = L;
    } // iter 

  if(L_) L_[0] = Lold;
  free(Sy);
}

