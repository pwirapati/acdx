/*
 * bremt.c - bootstrap resampling multiple testing
 *
 */

#include <stdlib.h>
#include <math.h>
#include "QSORT.h"

void
bremt(
  const int *dims_,
  double *theta,     // M test vs B bootstrap instances
  const double *T0_,
  double *SE,        // resmapling SE of theta
  double *PCER,
  double *FWER,
  double *FDR
  )
{
  const int M = dims_[0], B = dims_[1];
  double T0 = *T0_;

  int *nf = (int*) malloc( sizeof(int) * M * 3 );
  int *o = nf + M;
  int *ob = o + M;

  for(int i = 0; i < M; i++ )
    {
    SE[i] = 0;
    nf[i] = 0;
    PCER[i] = 1;
    }

  /*
   * scaling and centering for null distribution
   */

  // center by full theta and sum residual squares for SE
  for(int b = 1; b < B; b++ )
    {
    double *theta_b = theta + b*M;
    for(int i = 0; i < M; i++ )
      {
      if( !isfinite(theta_b[i]) ) continue;
      nf[i]++;
      theta_b[i] -= theta[i];
      SE[i] += theta_b[i]*theta_b[i];
      if(theta_b[i] > theta[i])
        PCER[i]++;
      }
    }
  for(int i = 0; i < M; i++ )
    {
    SE[i] = sqrt(SE[i]/nf[i]);
    PCER[i] /= nf[i]+1;
    }

  // turn theta to T by scaling;
  // convert-+/Inf/NA to -Inf to put last in sorting
  for(int b = 0; b < B; b++ )
    {
    double *theta_b = theta + b*M;
    for(int i = 0; i < M; i++ )
      {
      if(isfinite(theta_b[i]))
        theta_b[i] /= SE[i];
      else
        theta_b[i] = -INFINITY;
      }
    }

  /*
   *  step-down FWER and FDR
   *
   * Implement step-down maxT and ST-procedure
   * in Ge, Dudoit, Speed (2003) Test 12:1-77
   * algorithm in box 2 and 5, respectively.
   *
   */

  double *T = theta;
  
  for(int i = 0; i < M; i++ )
    {
    FWER[i] = 1;
    FDR[i] = 0;
    }

  for(int r = 0; r < M; r++ )
    o[r] = r;
  QSORT(int, M, int, o, T[_u] > T[_v]);

  double W0 = M, W0ave = 0;
  for(int r = 0; r < M; r++ )
    if( T[o[r]] <= T0 ) 
      { W0 -= r; break; }

  for(int b = 1; b < B; b++ )
    {
    double *Tb = T + b*M;
    for(int r = 0; r < M; r++ )
      ob[r] = r;
    QSORT( int, M, int, ob, Tb[_u] > Tb[_v] );
  
    // FWER
    double maxT = -INFINITY;
    for(int r = M-1; r >= 0; r-- )
      {
      if( Tb[o[r]] > maxT ) maxT = Tb[o[r]];
      if( maxT >= T[o[r]] )
        FWER[o[r]]++;
      }

    // FDR
    int q = 0;
    double V = 0;
    for(int r = 0;  r < M; r++ )
      {
      while( q < M-1 && Tb[ob[q]] >= T[o[r]] )
        { 
        V++;
        q++;
        }
      FDR[ o[r] ] += V;
      }

    double W0b = M;
    for(int r = 0; r < M; r++ )
      if( Tb[ob[r]] <= T0 ) 
        { W0b -= r; break; }
    W0ave += W0b;
    }

  FWER[o[0]] /= B;
  for(int r = 1; r < M; r++ )
    {
    FWER[ o[r] ] /= B;
    if( FWER[o[r-1]] > FWER[o[r]] )
      FWER[o[r]] = FWER[o[r-1]];
    }

  W0ave /= (B-1);
  FDR[o[M-1]] /= (B-1) * M;
  for(int r = M-2; r >= 0; r-- )
    {
    FDR[o[r]] /= (B-1) * (r+1) * W0ave/W0;
    if( FDR[o[r+1]] < FDR[o[r]])
      FDR[o[r]] = FDR[o[r+1]];
    }
  
  free(nf);
}
