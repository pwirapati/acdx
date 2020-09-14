/*
 * malm_acdx.c
 *
 *  - input: 
 *    - aggregate data from greg (sample x gene x ctype)
 *    - design matrix x
 *    - bootstrap resamples
 *
 *  - repeat for each bootstrap instances (first one is the original)
 *    - repeat for each cell-type:
 *      - normalize using malm11: alpha, beta
 *      - malm for each gene, with log(alpha) as offset [optional]
 *  
 *  - output
 *    - (alpha, beta, gamma, phi) x ctype x bootstraps
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "malm.h"

void
malm_acdx (
  const int *dim,     // dimensions
  uv *Y,              // n_sample x n_gene x n_ctype
  double *N,          // n_sample x n_ctype
  const double *Xt,   // n_coef x n_sample
  const int *Gi,      // n_disp
  const int *id_boot,    // n_sample: map sample to bootstrab block
  const int *B,       // n_bid x n_boot: resampled bid
  double *gamma,      // n_gene x n_ctype x n_boot
  double *phi0,       // n_gene x n_ctype x n_boot
  double *alpha,      // n_sample x n_ctype x n_boot
  double *beta,       // n_coef x n_gene x n_ctype x n_boot
  double *phi,        // n_disp x n_gene x n_ctype x n_boot
  const double *y0,
  const int *iopt11_,  
  const double *dopt11_,
  const int *iopt_,  
  const double *dopt_
  )
{
  int n_sample = dim[0];
  int n_gene = dim[1];
  int n_coef = dim[2];
  int n_disp = dim[3];   // dispersion groups
  int n_ctype = dim[4];
  int n_bid = dim[5];   // number of bootstrap blocks
  int n_boot = dim[6];   // number of bootstrap resamples

  // precompute sample weights for bootstrap, which is also
  // cell-type specific
  double *w = (double*)malloc(sizeof(double)
                  * n_sample * n_ctype * n_boot);
  for(int b = 0; b < n_boot; b++ )
    {
    const int *B_b = B + b * n_bid;
    double *w_b = w + b*n_sample*n_ctype;
    for(int k = 0; k < n_ctype; k++ )
      {
      const double *Nk = N + k*n_sample;
      double *w_bk = w_b + k*n_sample;
      for(int i = 0; i < n_sample; i++ )
        {
        w_bk[i] = 0;
        for(int a = 0; a < n_bid; a++ )
          if( id_boot[i] == B_b[a] ) w_bk[i]++;
        w_bk[i] *= ((Nk[i] > 1 ) ? 1 : 0);
        }
      }
    }

  
  for(int i = 0; i < n_sample * n_gene * n_ctype; i++ )
    Y[i].u += y0[0];
  
  const int n_thread = omp_get_max_threads();
  uv *E = (uv*)malloc( sizeof(uv) * n_sample * n_gene * n_thread );

  #pragma omp parallel for schedule(runtime)
  for(int b = 0; b < n_boot; b++ )
    {
    int t_id = omp_get_thread_num();
    const double *w_b = w + b * n_sample * n_ctype;
    double *alpha_b = alpha + b * n_sample * n_ctype;
    double *gamma_b = gamma + b * n_gene * n_ctype;
    double *phi0_b = phi0 + b * n_gene * n_ctype;
    double *beta_b = beta + b * n_coef * n_gene * n_ctype;
    double *phi_b = phi + b * n_disp * n_gene * n_ctype;
    uv *E_t = E + t_id * n_sample * n_gene;
    for(int k = 0; k < n_ctype; k++ )
      {
      uv *Yk = Y + k * n_sample * n_gene;
      const double *w_bk = w_b + k * n_sample;
      double *alpha_bk = alpha_b + k * n_sample;
      double *gamma_bk = gamma_b + k * n_gene;
      double *phi0_bk = phi0_b + k * n_gene;
      double *beta_bk = beta_b + k * n_coef * n_gene;
      double *phi_bk = phi_b + k * n_disp * n_gene;

      double L;
      malm11( dim, Yk, w_bk, E_t,
        alpha_bk, gamma_bk, phi0_bk, &L,
        iopt11_, dopt11_
        );
      
      malm(dim, Yk, w_bk, Xt, Gi, E_t,
        alpha_bk, beta_bk, phi_bk, NULL, iopt_, dopt_, NULL);
      } // ctype
    if(iopt_[0]==1) fputc('#',stderr);
    } // bootstrap

  if(iopt_[0]==1) fputc('\n',stderr);
  
  free(w);
  free(E);
  return;
}
