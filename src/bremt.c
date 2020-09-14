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
  const double *alpha_,
  double *SE,        // resmapling SE of theta
  double *PCER,
  double *FWER,
  double *FDP
  )
{
  const int M = dims_[0], B = dims_[1];
  double alpha = alpha_[0];
  if(alpha > 1 ) alpha = 1;
  if(alpha < 0 ) alpha = 0;

  int *nf = (int*) malloc( sizeof(int) * M * 2 );
  int *o = nf + M;
  double *maxT0 = (double*) malloc(sizeof(double) * B );

  for(int i = 0; i < M; i++ )
    {
    SE[i] = 0;
    nf[i] = 0;
    o[i] = i;
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

  // turn theta to T by scaling; -+/Inf/NA to -Inf to put last in sorting
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
   *  FWER
   */

  // sort to obtain maxT distribution (on the first row)
  for(int b = 1; b < B; b++ )
    {
    double *theta_b = theta + b*M;
    QSORT( int, M, double, theta_b, _u > _v );
    }

  // copy and sort the maximum
  int Bf = 0;
  for(int b = 1; b < B; b++ )
    {
    double T = theta[b*M];
    if(isfinite(T)) Bf++;
    maxT0[b] = T;
    }
  QSORT(int, Bf, double, maxT0, _u > _v );

  for(int i = 0; i < M; i++ )
    {
    FWER[i] = 1;
    for(int b = 0; b < Bf; b++ )
      if( theta[i] < maxT0[b] ) 
        FWER[i]++;
    FWER[i] /= Bf+1;
    }  


  /*
   *  FDP
   */
  // obtain empircal quantile of T at p = 1-1/M,1-2/M,1-3/M, ...
  
  // put the values on the second column
  double *Tnull = theta + M;
  Tnull[0] = maxT0[ (int)(alpha * (Bf-1))];  // first one already available
  for(int i = 1; i < M; i++)
    {
    Bf = 0;
    for(int b = 1; b < B; b++ )
      {
      double T = theta[i + b*M];
      if(isfinite(T)) Bf++;
      maxT0[b] = T;
      }
    if( Bf == 0 ) 
      {
      Tnull[i] = Tnull[i-1]; // conservative truncation
      continue;
      }
    QSORT(int, Bf, double, maxT0, _u > _v );
    Tnull[i] = maxT0[ (int)(alpha * (Bf-1))];
    }

  // ordering of the observed T (for the number of 'discovery')
  QSORT( int, M, int, o, theta[_u] > theta[_v] );
  for(int i = 0; i < M; i++ )
    {
    double T = theta[o[i]];
    int j;
    for(j = 0; j < M && Tnull[j] > T; j++ )
      j++;
    FDP[o[i]] = (double)j/(i+1);
    }
  for(int i = M-2; i >= 0; i-- ) // step-up
    if( FDP[o[i]] > FDP[o[i+1]] ) FDP[o[i]] = FDP[o[i+1]];

  free(maxT0);
  free(nf);
}
