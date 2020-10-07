/*
 * avediff.c - averaged-coefficient difference
 *
 */

#include <math.h>

void
avediff(
  const int *dims_,
  const double *beta,
  int *nH, int *Hctype, int *Hcoef,
  int *nL, int *Lctype, int *Lcoef,
  double *theta,
  int *na_rm
  )
{
  const int n_coef=dims_[0];
  const int n_gene=dims_[1];
  const int n_ctype=dims_[2];
  const int n_boot=dims_[3];

  for(int b = 0; b < n_boot; b++ )
    {
    double *theta_b = theta + b*n_gene;
    const double *beta_b = beta + b*n_coef*n_gene*n_ctype;
    for(int j = 0; j < n_gene; j++ )
      {
      const double *beta_bj = beta_b + j*n_coef;

      double aveH = 0, n = 0;
      for(int h = 0; h < nH[0]; h++ )
        {
        const double beta_h = beta_bj[ Hctype[h] * n_coef * n_gene + Hcoef[h]]; 
        if( !na_rm && isnan(beta_h))
          { aveH = NAN; break; }
        aveH += beta_h;
        n++;
        }
      aveH /= n;

      double aveL = 0;
      n = 0;
      for(int l = 0; l < nL[0]; l++ )
        {
        const double beta_l = beta_bj[ Lctype[l] * n_coef * n_gene + Lcoef[l]];
        if( !na_rm && isnan(beta_l))
          { aveL = NAN; break; }
        aveL += beta_l;
        n++;
        }
      aveL /= n;
      
      theta_b[j] = aveH - aveL;
      }
    }
}
