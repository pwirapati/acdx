/*
 * nglm.c 
 * 
 * - response variable are pairs of mean and variance
 * - multiple responses
 *    - main design 'X` replicated to all variables
 *    - linked by intercept for a group of observations
 * - quadratic variance-mean model (any component can be selectively disabled)
 *   
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "nglm.h"
#include "nchol.h"

// init values for gamma model
// - fit normal model (ignore Y.v) in shifted log space
int
nglm_init_gamma_ss (
  const int ni,
  const int nj,
  const int np,
  const int nq,

  const double *w,
  const uv *Y,
  const int *Xp,
  const int *Zq,

  double *alpha,          // 1*nq
  double *v_alpha,        // 1*1*nq    [Not used]
  double *beta,           // np*nj
  double *v_beta,         // np*np*nj  [Not used]
  double *phi,            // 1*nj
  uv *E,                  // ni*nj
  
  const double shift,     // response = log( y + shift )

  const double tol,
  const int itermax       // one is enough 
  )
{
  const int nthread = omp_get_max_threads();
  const int szws = 2*nq + np + np;
  double *ws = (double*)malloc(sizeof(double)*nthread* szws);

  double sumw = 0;
  for(int i = 0; i < ni; i++ )
    sumw += w[i];

  // initialize all params to zero
  //
  for(int ij = 0; ij < ni*nj; ij++ )
    E[ij].u = log( Y[ij].u + shift );
  for(int q = 0; q < nq; q++ )
    alpha[q] = 0;
  for(int pj = 0; pj < np*nj; pj++ )
    beta[pj] = 0;

  // iterate between alpha and beta updates
  for(int iter = 0; iter < itermax; iter++ )
    {
    double norm1_D_eta = 0;

    /*
     *  update alpha
     */
    for(int tid = 0; tid < nthread; tid++ )
      {
      double *SHz = ws + tid*szws;
      for(int q = 0; q < 2*nq; q++ )
        SHz[q] = 0;
      }

    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sz = ws + tid*szws;
      double *Hz = Sz + nq;
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        Sz[ Zq[i] ] += w[i]*E_j[i].u;
        Hz[ Zq[i] ] += w[i];
        }
      }

    double *Sz0 = ws;
    double *Hz0 = Sz0 + nq;
    for(int tid = 1; tid < nthread; tid++ )
      {
      double *Szt = ws + tid*szws;
      double *Hzt = Szt + nq;
      for(int q = 0; q < nq; q++ )
        {
        Sz0[q] += Szt[q];
        Hz0[q] += Hzt[q];
        }
      }

    double mean_Dalpha = 0;
    for(int q = 0; q < nq; q++ )
      {
      if( Hz0[q] == 0 ) Sz0[q] = 0;
      else Sz0[q] /= Hz0[q];     // now it is the change of alpha
      mean_Dalpha += Sz0[q];
      }
    mean_Dalpha /= nq;
    for(int q = 0; q < nq; q++ ) 
      {
      Sz0[q] -= mean_Dalpha;    // sum-to-zero constrain
      alpha[q] += Sz0[q];       // update alpha
      }

    #pragma omp parallel for reduction(+: norm1_D_eta) schedule(runtime)
    for(int j = 0; j < nj; j++ ) 
      {
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        E_j[i].u -= Sz0[Zq[i]];
        norm1_D_eta += w[i]*fabs(Sz0[Zq[i]]);
        }
      }

    /*
     *  update beta
     */
    #pragma omp parallel for reduction(+: norm1_D_eta) schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sx = ws + tid*szws + 2*nq;
      double *Hxx = Sx + np;

      for(int p = 0; p < np; p++ )
        Sx[p] = Hxx[p] = 0;
      
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        Sx[Xp[i]] += w[i]*E_j[i].u;
        Hxx[Xp[i]] += w[i];
        }

      for(int p = 0; p < np; p++ )
        Sx[p] /= Hxx[p];

      double *beta_j = beta + j*np;
      for(int p = 0; p < np; p++ )
        beta_j[p] += Sx[p];
      for(int i = 0; i < ni; i++ )
        {
        E_j[i].u -= Sx[Xp[i]];
        norm1_D_eta += w[i]*fabs(Sx[Xp[i]]);
        }
      }

    //fprintf(stderr,"[%4d]: %g\n", iter,norm1_D_eta/(sumw*nj));
    if( norm1_D_eta/(sumw*nj) < tol ) break;
    }

  // phi per gene estimated as variance in log scale
  // replace E[ij].u by mu[ij] = exp(eta[ij])
  // replace E[ij].v by phi[j]*mu[ij]*mu[ij]
  const double df = sumw - np - (double)nq/nj; 
  #pragma omp parallel for schedule(runtime)
  for(int j = 0; j < nj; j++ )
    {
    phi[j] = 0;
    uv *E_j = E + j*ni;
    const uv *Y_j = Y + j*ni;
    for(int i = 0; i < ni; i++ )
      phi[j] += w[i]*E_j[i].u*E_j[i].u;
    phi[j] /= df;

    for(int i = 0; i < ni; i++ )
      {
      double eta_ij = log( Y_j[i].u + shift ) - E_j[i].u;
      E_j[i].u = exp(eta_ij);
      E_j[i].v = Y_j[i].v + phi[j]*E_j[i].u*E_j[i].u; // total var added
      }
    }

  free(ws);
  return 0;
}
// init values for gamma model
// - fit normal model (ignore Y.v) in shifted log space
int
nglm_init_gamma_sd (
  const int ni,
  const int nj,
  const int np,
  const int nq,

  const double *w,
  const uv *Y,
  const double *Xt,
  const int *Zq,

  double *alpha,          // 1*nq
  double *v_alpha,        // 1*1*nq    [Not used]
  double *beta,           // np*nj
  double *v_beta,         // np*np*nj  [Not used]
  double *phi,            // 1*nj
  uv *E,                  // ni*nj
  
  const double shift,     // response = log( y + shift )

  const double tol,
  const int itermax       // one is enough 
  )
{
  const int nh = np*(np+1)/2;
  const int nthread = omp_get_max_threads();
  const int szws = 2*nq + np + nh;
  double *ws = (double*)malloc(sizeof(double)*nthread* szws);

  double sumw = 0;
  for(int i = 0; i < ni; i++ )
    sumw += w[i];

  // initialize all params to zero
  //
  for(int ij = 0; ij < ni*nj; ij++ )
    E[ij].u = log( Y[ij].u + shift );
  for(int q = 0; q < nq; q++ )
    alpha[q] = 0;
  for(int pj = 0; pj < np*nj; pj++ )
    beta[pj] = 0;

  // iterate between alpha and beta updates
  for(int iter = 0; iter < itermax; iter++ )
    {
    double norm1_D_eta = 0;

    /*
     *  update alpha
     */
    for(int tid = 0; tid < nthread; tid++ )
      {
      double *SHz = ws + tid*szws;
      for(int q = 0; q < 2*nq; q++ )
        SHz[q] = 0;
      }

    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sz = ws + tid*szws;
      double *Hz = Sz + nq;
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        Sz[ Zq[i] ] += w[i]*E_j[i].u;
        Hz[ Zq[i] ] += w[i];
        }
      }

    double *Sz0 = ws;
    double *Hz0 = Sz0 + nq;
    for(int tid = 1; tid < nthread; tid++ )
      {
      double *Szt = ws + tid*szws;
      double *Hzt = Szt + nq;
      for(int q = 0; q < nq; q++ )
        {
        Sz0[q] += Szt[q];
        Hz0[q] += Hzt[q];
        }
      }

    double mean_Dalpha = 0;
    for(int q = 0; q < nq; q++ )
      {
      if( Hz0[q] == 0 ) Sz0[q] = 0;
      else Sz0[q] /= Hz0[q];     // now it is the change of alpha
      mean_Dalpha += Sz0[q];
      }
    mean_Dalpha /= nq;
    for(int q = 0; q < nq; q++ ) 
      {
      Sz0[q] -= mean_Dalpha;    // sum-to-zero constrain
      alpha[q] += Sz0[q];       // update alpha
      }

    #pragma omp parallel for reduction(+: norm1_D_eta) schedule(runtime)
    for(int j = 0; j < nj; j++ ) 
      {
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        E_j[i].u -= Sz0[Zq[i]];
        norm1_D_eta += w[i]*fabs(Sz0[Zq[i]]);
        }
      }

    /*
     *  update beta
     */
    #pragma omp parallel for reduction(+: norm1_D_eta) schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sx = ws + tid*szws + 2*nq;
      double *Hx = Sx + np;

      for(int p = 0; p < np; p++ )
        Sx[p] = 0;
      for(int h = 0; h < nh; h++ )
        Hx[h] = 0;
      
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        const double *X_i = Xt + i*np;
        int h = 0;
        for(int p = 0; p < np; p++ )
          {
          Sx[p] += w[i] * X_i[p] * E_j[i].u;
          for(int pp = 0; pp <= p; pp++, h++ )
            Hx[h] += w[i] * X_i[p] * X_i[pp];
          }
        }

      nchol( np, Hx, NCHOL_AUTOTOL );
      nchol_fsub( np, Hx, Sx );
      nchol_bsub( np, Hx, Sx );
      
      double *beta_j = beta + j*np;
      for(int p = 0; p < np; p++ )
        beta_j[p] += Sx[p];
      for(int i = 0; i < ni; i++ )
        {
        const double *X_i = Xt + i*np;
        double D_eta = 0;
        for(int p = 0; p < np; p++ )
          D_eta -= Sx[p] * X_i[p];
        E_j[i].u += D_eta;
        norm1_D_eta += w[i]*fabs(D_eta);
        }
      }

    fprintf(stderr,"[%4d]: %g\n", iter,norm1_D_eta/(sumw*nj));
    if( norm1_D_eta/(sumw*nj) < tol ) break;
    }

  // phi per gene estimated as variance in log scale
  // replace E[ij].u by mu[ij] = exp(eta[ij])
  // replace E[ij].v by phi[j]*mu[ij]*mu[ij]
  const double df = sumw - np - (double)nq/nj; 
  #pragma omp parallel for schedule(runtime)
  for(int j = 0; j < nj; j++ )
    {
    phi[j] = 0;
    uv *E_j = E + j*ni;
    const uv *Y_j = Y + j*ni;
    for(int i = 0; i < ni; i++ )
      phi[j] += w[i] * E_j[i].u*E_j[i].u;
    phi[j] /= df;

    for(int i = 0; i < ni; i++ )
      {
      double eta_ij = log( Y_j[i].u + shift ) - E_j[i].u;
      E_j[i].u = exp(eta_ij);
      E_j[i].v = Y_j[i].v + phi[j]*E_j[i].u*E_j[i].u; // total var added
      }
    }

  free(ws);
  return 0;
}

//
// nglm_gamma with sparse Z and sparse X
//

int
nglm_gamma_ss (
  const int ni,
  const int nj,
  const int np,
  const int nq,

  const double *w,
  const uv *Y,
  const int *Xp,          // ni
  const int *Zq,          // ni

  double *alpha,          // 1*nq
  double *v_alpha,        // 1*1*nq
  double *beta,           // np*nj
  double *v_beta,         // np*nj
  double *phi,            // 1*nj
  uv *E,                  // ni*nj
  
  const double shift,     // response = log( y + shift )

  const double tol,
  const int itermax
  )
{
  // initialize
  nglm_init_gamma_ss(
    ni, nj, np, nq,
    w, Y, Xp, Zq,
    alpha, NULL, beta, NULL, phi, E,
    0.5,
    tol, 1
    );

  const int nthread = omp_get_max_threads();
  const int szws = 2*nq + np + np;
  double *ws = (double*)malloc(sizeof(double)*nthread* szws);

  double sumw = 0;
  for(int i = 0; i < ni; i++ )
    sumw += w[i];

  // iterate between alpha and beta updates
  for(int iter = 0; iter < itermax; iter++ )
    {
    double norm1_D_mu = 0;

    /*
     *  update alpha
     */
    for(int tid = 0; tid < nthread; tid++ )
      {
      double *SHz = ws + tid*szws;
      for(int q = 0; q < 2*nq; q++ )
        SHz[q] = 0;
      }

    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sz = ws + tid*szws;
      double *Hz = Sz + nq;
      uv *E_j = E + j*ni;
      const uv *Y_j = Y + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        if(w[i]==0) continue;
        double W = w[i] * E_j[i].u/E_j[i].v;
        Sz[ Zq[i] ] += W*(Y_j[i].u - E_j[i].u);
        Hz[ Zq[i] ] += W*E_j[i].u;
        }
      }

    double *Sz0 = ws;
    double *Hz0 = Sz0 + nq;
    for(int tid = 1; tid < nthread; tid++ )
      {
      double *Szt = ws + tid*szws;
      double *Hzt = Szt + nq;
      for(int q = 0; q < nq; q++ )
        {
        Sz0[q] += Szt[q];
        Hz0[q] += Hzt[q];
        }
      }

    double mean_Dalpha = 0;
    for(int q = 0; q < nq; q++ )
      {
      v_alpha[q] = 1/Hz0[q];
      if(Hz0[q]==0) Sz0[q] = 0;
      else Sz0[q] /= Hz0[q];       // the change of alpha
      mean_Dalpha += Sz0[q];
      }
    mean_Dalpha /= nq;
    for(int q = 0; q < nq; q++ ) 
      {
      Sz0[q] -= mean_Dalpha;    // sum-to-zero constrain
      alpha[q] += Sz0[q];       // update alpha
      Sz0[q] = exp(Sz0[q]);     // for modifying mu
      }

    // update mu
    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < nj; j++ ) 
      {
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        E_j[i].u *= Sz0[Zq[i]];
      }

    /*
     *  update beta
     */
    #pragma omp parallel for reduction(+:norm1_D_mu) schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sx = ws + tid*szws + 2*nq;
      double *Hxx = Sx + np;

      for(int p = 0; p < np; p++ )
        Sx[p] = 0;
      for(int p = 0; p < np; p++ )
        Hxx[p] = 0;
      
      uv *E_j = E + j*ni;
      const uv *Y_j = Y + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        if(w[i]==0) continue;
        double Wx = w[i]*E_j[i].u/E_j[i].v;
        const int pi = Xp[i];
        Sx[pi] += Wx * (Y_j[i].u - E_j[i].u);
        Wx *= E_j[i].u;
        Hxx[pi] += Wx;
        }

      for(int p = 0; p < np; p++ )
        if( Hxx[p] == 0 ) Sx[p] = 0;
        else Sx[p] /= Hxx[p];

      double *beta_j = beta + j*np;
      for(int p = 0; p < np; p++ )
        beta_j[p] += Sx[p];
      for(int i = 0; i < ni; i++ )
        {
        double eta = alpha[Zq[i]] + beta_j[Xp[i]];
        double mu = exp(eta);
        norm1_D_mu += w[i]*fabs(E_j[i].u - mu);
        E_j[i].u = mu;
        }

      /*
       * update phi
       */
      double Hg = 0, Sg = 0;
      double varbiascorrect = (double)sumw/(sumw-np-nq/nj);
      for(int i = 0; i < ni; i++ )
        {
        if(w[i]==0) continue;
        double W = w[i] * E_j[i].u * E_j[i].u / E_j[i].v;
        double D = Y_j[i].u - E_j[i].u; 
        double D2 = D*D * varbiascorrect;
        Sg += W * D2 / E_j[i].v;
        Hg += W;
        }
      phi[j] *= Sg/Hg;
      for(int i = 0; i < ni; i++ )
        E_j[i].v = Y_j[i].v + phi[j]*E_j[i].u*E_j[i].u;
      }

    //fprintf(stderr,"[%4d]: %g\n", iter,norm1_D_mu/(sumw*nj));
    fputc('*',stderr);
    if( norm1_D_mu/(sumw*nj) < tol ) break;
    }
  fputc('\n',stderr);

  // vcov  beta
  #pragma omp parallel for schedule(runtime)
  for(int j = 0; j < nj; j++ )
    {
    const int tid = omp_get_thread_num();
    double *Hxx = ws + tid*szws;
    for(int p = 0; p < np; p++ )
      Hxx[p] = 0;
    const uv* E_j = E + j*ni;
    for(int i = 0; i < ni; i++ )
      {
      if(w[i]==0) continue;
      Hxx[Xp[i]] += w[i] * E_j[i].u*E_j[i].u/E_j[i].v;
      }

    double *v_beta_j = v_beta + j*np;
    for(int p = 0; p < np; p++ )
      v_beta_j[p] = 1/Hxx[p];

    }

  free(ws);
  return 0;
}


//
// nglm_gamma with sparse Z and dense X
//

int
nglm_gamma_sd (
  const int ni,
  const int nj,
  const int np,
  const int nq,

  const double *w,
  const uv *Y,
  const double *Xt,
  const int *Zq,

  double *alpha,          // 1*nq
  double *v_alpha,        // 1*1*nq
  double *beta,           // np*nj
  double *v_beta,         // np*np*nj
  double *phi,            // 1*nj
  uv *E,                  // ni*nj
  
  const double shift,     // response = log( y + shift )

  const double tol,
  const int itermax
  )
{
  // initialize 
  nglm_init_gamma_sd(
    ni, nj, np, nq,
    w, Y, Xt, Zq,
    alpha, NULL, beta, NULL, phi, E,
    0.5,
    tol, 1
    );

  const int nh = np*(np+1)/2;
  const int nthread = omp_get_max_threads();
  const int szws = 2*nq + np + nh;
  double *ws = (double*)malloc(sizeof(double)*nthread* szws);

  double sumw = 0;
  for(int i = 0; i < ni; i++ )
    sumw += w[i];

  // iterate between alpha and beta updates
  for(int iter = 0; iter < itermax; iter++ )
    {
    double norm1_D_mu = 0;

    /*
     *  update alpha
     */
    for(int tid = 0; tid < nthread; tid++ )
      {
      double *SHz = ws + tid*szws;
      for(int q = 0; q < 2*nq; q++ )
        SHz[q] = 0;
      }

    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sz = ws + tid*szws;
      double *Hz = Sz + nq;
      uv *E_j = E + j*ni;
      const uv *Y_j = Y + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        if(w[i] == 0 ) continue;
        double W = w[i]*E_j[i].u/E_j[i].v;
        Sz[ Zq[i] ] += W*(Y_j[i].u - E_j[i].u);
        Hz[ Zq[i] ] += W*E_j[i].u;
        }
      }

    double *Sz0 = ws;
    double *Hz0 = Sz0 + nq;
    for(int tid = 1; tid < nthread; tid++ )
      {
      double *Szt = ws + tid*szws;
      double *Hzt = Szt + nq;
      for(int q = 0; q < nq; q++ )
        {
        Sz0[q] += Szt[q];
        Hz0[q] += Hzt[q];
        }
      }

    double mean_Dalpha = 0;
    for(int q = 0; q < nq; q++ )
      {
      v_alpha[q] = 1/Hz0[q];
      if(Hz0[q]==0) Sz0[q] = 0;
      else Sz0[q] /= Hz0[q];       // the change of alpha
      mean_Dalpha += Sz0[q];
      }
    mean_Dalpha /= nq;
    for(int q = 0; q < nq; q++ ) 
      {
      Sz0[q] -= mean_Dalpha;    // sum-to-zero constrain
      alpha[q] += Sz0[q];       // update alpha
      Sz0[q] = exp(Sz0[q]);     // for modifying mu
      }

    // update mu
    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < nj; j++ ) 
      {
      uv *E_j = E + j*ni;
      for(int i = 0; i < ni; i++ )
        E_j[i].u *= Sz0[Zq[i]];
      }

    /*
     *  update beta
     */
    #pragma omp parallel for reduction(+:norm1_D_mu) schedule(runtime)
    for(int j = 0; j < nj; j++ )
      {
      const int tid = omp_get_thread_num();
      double *Sx = ws + tid*szws + 2*nq;
      double *Hx = Sx + np;

      for(int p = 0; p < np; p++ )
        Sx[p] = 0;
      for(int h = 0; h < nh; h++ )
        Hx[h] = 0;
      
      uv *E_j = E + j*ni;
      const uv *Y_j = Y + j*ni;
      for(int i = 0; i < ni; i++ )
        {
        if(w[i]==0) continue;
        const double *X_i = Xt + i*np;
        int h = 0;
        for(int p = 0; p < np; p++ )
          {
          double wx = w[i]*E_j[i].u/E_j[i].v * X_i[p];
          Sx[p] += wx * (Y_j[i].u - E_j[i].u);
          wx *= E_j[i].u;
          for(int pp = 0; pp <= p; pp++, h++ )
            Hx[h] += wx * X_i[pp];
          }
        }

      nchol( np, Hx, NCHOL_AUTOTOL );
      nchol_fsub( np, Hx, Sx );
      nchol_bsub( np, Hx, Sx );
      
      double *beta_j = beta + j*np;
      for(int p = 0; p < np; p++ )
        beta_j[p] += Sx[p];
      for(int i = 0; i < ni; i++ )
        {
        const double *X_i = Xt + i*np;
        double eta = alpha[Zq[i]];
        for(int p = 0; p < np; p++ )
          eta += beta_j[p] * X_i[p];
        double mu = exp(eta);
        norm1_D_mu += w[i]*fabs(E_j[i].u - mu);
        E_j[i].u = mu;
        }

      /*
       * update phi
       */
      double Hg = 0, Sg = 0;
      double varbiascorrect = sumw/(sumw-np-nq/nj);
      for(int i = 0; i < ni; i++ )
        {
        if( w[i] == 0 ) continue;
        double W = w[i]*E_j[i].u * E_j[i].u / E_j[i].v;
        double D = Y_j[i].u - E_j[i].u; 
        double D2 = D*D * varbiascorrect;
        Sg += W * D2 / E_j[i].v;
        Hg += W;
        }
      phi[j] *= Sg/Hg;
      for(int i = 0; i < ni; i++ )
        E_j[i].v = Y_j[i].v + phi[j]*E_j[i].u*E_j[i].u;
      }

    fprintf(stderr,"[%4d]: %g\n", iter,norm1_D_mu/(sumw*nj));
    if( norm1_D_mu/(sumw*nj) < tol ) break;
    }

  // vcov  beta
  #pragma omp parallel for schedule(runtime)
  for(int j = 0; j < nj; j++ )
    {
    const int tid = omp_get_thread_num();
    double *Hx = ws + tid*szws;
    for(int h = 0; h < nh; h++ )
      Hx[h] = 0;
    const uv* E_j = E + j*ni;
    for(int i = 0; i < ni; i++ )
      {
      if( w[i] == 0 ) continue;
      const double *X_i = Xt + i*np;
      int h = 0;
      for(int p = 0; p < np; p++ )
        {
        double Wx = w[i] * E_j[i].u*E_j[i].u/E_j[i].v *  X_i[p];
        for(int pp = 0; pp <= p; pp++, h++ )
          Hx[h] += Wx * X_i[pp];
        }
      }

    nchol( np, Hx, NCHOL_AUTOTOL );

    double *v_beta_j = v_beta + j*np*np;
    for(int p = 0; p < np; p++ )
      {
      double *v_beta_jp =  v_beta_j + np*p;
      for(int pp = 0; pp < np; pp++ )
        v_beta_jp[pp] = (pp == p ? 1: 0);
      nchol_fsub( np, Hx, v_beta_jp );
      nchol_bsub( np, Hx, v_beta_jp );
      }
    }

  free(ws);
  return 0;
}
