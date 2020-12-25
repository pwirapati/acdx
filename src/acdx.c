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
 *    - (alpha, beta, beta0, phi) x ctype x bootstraps
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
  double *beta0,      // n_gene x n_ctype x n_boot
  double *phi0,       // n_gene x n_ctype x n_boot
  double *alpha,      // n_sample x n_ctype x n_boot
  double *beta,       // n_coef x n_gene x n_ctype x n_boot
  double *phi,        // n_disp x n_gene x n_ctype x n_boot
  const int *norm_method,
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


  // allocate thread-specific storage
  //
  const int n_thread = omp_get_max_threads();
  uv *E = (uv*)malloc( sizeof(uv) * n_sample * n_gene * n_thread );

  uv *Yg, *Eg;
  if( *norm_method == 2 )
    {
    Eg = (uv*)malloc( sizeof(uv) * n_sample * n_gene * n_ctype * n_thread);
    Yg = (uv*)malloc( sizeof(uv) * n_sample * n_gene * n_ctype );
    for(int k = 0; k < n_ctype; k++ )
      for(int j = 0; j < n_gene; j++ )
        for(int i = 0; i < n_sample; i++ )
          Yg[i + k*n_sample + j*n_sample*n_ctype ] =
             Y[i + j*n_sample + k*n_sample*n_gene ];
    }

  int n_done = 0, chunk=ceil(n_boot/100.0);

  #pragma omp parallel for schedule(runtime)
 for(int b = 0; b < n_boot; b++ )
    {
    int t_id = omp_get_thread_num();
    const double *w_b = w + b * n_sample * n_ctype;
    double *alpha_b = alpha + b * n_sample * n_ctype;
    double *beta0_b = beta0 + b * n_gene * n_ctype;
    double *phi0_b = phi0 + b * n_gene * n_ctype;
    double *beta_b = beta + b * n_coef * n_gene * n_ctype;
    double *phi_b = phi + b * n_disp * n_gene * n_ctype;
    uv *E_t = E + t_id * n_sample * n_gene;

    if( *norm_method == 2 ) // global
      {
      int dimg[2] = { n_sample * n_ctype, n_gene };
      malm11( dimg, Yg, w_b, Eg + t_id * n_sample * n_ctype * n_gene,
        alpha_b, beta0_b, phi0_b, NULL,
        iopt11_,dopt11_ );
      }
    else if( *norm_method == 1 ) // per cell type
      {
      for(int k = 0; k < n_ctype; k++ )
        {
        uv *Yk = Y + k * n_sample * n_gene;
        const double *w_bk = w_b + k * n_sample;
        double *alpha_bk = alpha_b + k * n_sample;
        double *beta0_bk = beta0_b + k * n_gene;
        double *phi0_bk = phi0_b + k * n_gene;

        malm11( dim, Yk, w_bk, E_t,
          alpha_bk, beta0_bk, phi0_bk, NULL,
          iopt11_, dopt11_
          );
        }
      }
    else
      ;  // no normalization

    if(n_coef > 0 )
      for(int k = 0; k < n_ctype; k++ )
        {
        uv *Yk = Y + k * n_sample * n_gene;
        const double *w_bk = w_b + k * n_sample;
        double *alpha_bk = alpha_b + k * n_sample;

        double *beta_bk = beta_b + k * n_coef * n_gene;
        double *phi_bk = phi_b + k * n_disp * n_gene;

        malm(dim, Yk, w_bk, Xt, Gi, E_t,
          alpha_bk, beta_bk, phi_bk, NULL,
          iopt_, dopt_, NULL);

        }

    #pragma omp critical
      {
      n_done++;
      if(iopt_[1]==1 && (n_done % chunk == 0))
        fprintf(stderr,"\rbootstrap completed: %d (%.1f%%)",
          n_done, (double)(n_done)*100/n_boot);
      }
    } // bootstrap

  if(iopt_[1]==1) fputc('\n',stderr);

  if( *norm_method == 2 )
    {
    for(int k = 1; k < n_ctype; k++ )
      for(int j = 0; j < n_gene; j++ )
        beta0[j + k*n_gene] = beta0[j];
    free(Yg);
    free(Eg);
    }
  free(w);
  free(E);
  return;
}
