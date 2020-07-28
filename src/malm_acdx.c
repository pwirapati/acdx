/*
 * malm_acdx.c
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <omp.h>

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

// -2 * normal ("pseudo"-)log-likelihood + C
//
// Data points with infinite expected variance is ignored (assumed
// to be a constant)
//
void
pll ( 
  const int *nrow_,
  const int *ncol_,
  const double *w,
  const uv *Y,
  const uv *E,
  double *Qd_, double *Qe_
  )
{
  const int n = *nrow_, m = *ncol_;
  double Qd = 0, Qe = 0;
  #pragma omp parallel for reduction(+:Qd,Qe) schedule(runtime) 
  for(int j = 0; j < m; j++ )
    {
    const uv* Y_j = Y + j * n;
    const uv* E_j = E + j * n;

    for(int i = 0; i < n; i++ )
      {
      if(w[i] == 0 ) continue;
      double d = Y_j[i].u - E_j[i].u;
      Qd += w[i]*d*d/E_j[i].v;
      if(isfinite(E_j[i].v)) Qe += w[i]*log(E_j[i].v);
      }
    }
  *Qd_ = Qd;
  *Qe_ = Qe;
}

void
malm_acdx (
  const int *dims,

  const double *w,  // [n_agg]
  const uv *Y,      // (y,v) x n_ac x n_gene

  // mean submodel
  const int *X1,      // intercepts
  const int *Xf,      // factors (sum-to-zero constrained)
  const double *Xc,   // continuous covariates
  
  const int *Zf,      // per sample scaling (sum-to-zero constrained)
  
  // dispersion submodel
  const int *G1,    // intercepts

  // parameters
  double *beta,    // 
  double *alpha,   // 
  double *phi,     //

  // predicted values
  uv *E,           // (mu,sigma2) x n_ac x n_gene

  const double *tol_,
  const int *limits_,
  const int *options,
  int *status_
  )
{
  const int n_ac = dims[0];
  const int n_gene = dims[1];
  const int n_Zf = dims[2]; // library scaling factor
  const int n_X1 = dims[3]; // intercept factors (0 or 1)
  const int n_Xf = dims[4]; // sum-to-zero factors
  const int n_Xc = dims[5]; // continuous covariates
  const int n_G1 = dims[6];   // dispersion intercepts (0 or 1)
  const double iter_tol = tol_[0];
  const double bt_tol = tol_[1];
  const double inv_tol = tol_[2];
  const int iter_max = limits_[0];
  const int bt_max = limits_[1];
  const int wssize_min = limits_[2]; // workspace minimum size (# of doubles)
  const int default_init = options[0];
  const int verbose = options[1];


  // the number of parameters per gene
  int k_Xf_[n_Xf], k_Xf = 0;
  for(int p = 0; p < n_Xf; p++ )
    {
    k_Xf_[p] = nlevels(n_ac, Xf + n_ac * p);
    k_Xf += k_Xf_[p];
    };

  const int k_X1 = nlevels(n_ac, X1 );
  const int k_Zf = nlevels(n_ac, Zf );
  const int k_G1 = nlevels(n_ac, G1 );

  const int sz_beta = k_X1 + k_Xf + n_Xc;

  int szws = 2*k_Zf;
  if( szws < 2*k_X1 ) szws = 2*k_X1;
  if( szws < 2*k_Xf ) szws = 2*k_Xf;
  if( szws < 2*n_Xc ) szws = 2*n_Xc;
  if( szws < 2*k_G1 ) szws = 2*k_G1;
  if( szws < wssize_min ) szws = wssize_min;

  int n_thread = omp_get_max_threads();
  double* ws[n_thread];
  for(int t = 0; t < n_thread; t++ )
    {
    ws[t] = (double*)malloc(sizeof(double)*szws);
    if(!ws[t]) 
      { status_[0] = -1; return; }
    }

  /*
   *  initialize parameters and expectations
   *
   */
  if( default_init && n_X1 > 0 )
    {
    const double min_phi = 0.01;
    
    // library scale factor: init to all zeroes
    dset( k_Zf, alpha, 0 );

    // gene-specific intercepts
    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < n_gene; j++ )
      {
      const uv* Y_j = Y + j*n_ac;
      uv* E_j = E + j*n_ac;
      double* beta_j = beta + j*sz_beta;
      double* phi_j = phi + j*k_G1;

      // mean part

      double *A = ws[omp_get_thread_num()];
      double *B = A + k_X1;
      dset( 2*k_X1, A, 0 );

      for(int i = 0; i < n_ac; i++ )
        {
        if(w[i] == 0 ) continue;
        A[X1[i]] += w[i]/Y_j[i].v * Y_j[i].u; 
        B[X1[i]] += w[i]/Y_j[i].v;
        }
      for(int l = 0; l < k_X1; l++ )
        {
        if( B[l] > 0 ) A[l] /= B[l];
        else A[l] = 0;
        if(B[l] > 0) beta_j[l] = log(A[l]);
        }
      for(int i = 0; i < n_ac; i++ )
        {
        if(w[i]==0) continue;
        E_j[i].u = A[X1[i]];
        }

      // init all others to zero
      dset( sz_beta, beta + k_X1, 0 );

      // dispersion part
      if( n_G1 > 0 )
        {
        B = A + k_G1;
        dset( 2*k_G1, A, 0 );
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i]==0) continue;
          double d = Y_j[i].u - E_j[i].u;
          A[G1[i]] += w[i] * d * d;
          B[G1[i]] += w[i] * Y_j[i].v;
          }
        for(int l = 0; l < k_G1; l++ )
          {
          phi_j[l] = A[l]/B[l];
          if(phi_j[l] < min_phi )
            phi_j[l] = min_phi;
          }
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i] == 0 ) continue;
          E_j[i].v = Y_j[i].v + phi_j[G1[i]] * E_j[i].u * E_j[i].u;
          }
        }
      else  // FEMA (use only cell-level variance)
        {
        for(int i = 0; i < n_ac; i++ )
          E_j[i].v = Y_j[i].v;
        }
      }
    }
  else 
    mess("X intercept not specified, user init assumed\n"); 

  /*
   *  Model-Fitting Iteration
   */

  double Qd = 0, Qe = 0;
  pll( &n_ac, &n_gene, w, Y, E, &Qd, &Qe );
  double Qold = Qd + Qe;

  int nbt = 0;
  int iter;
  for(iter = 1; iter <= iter_max; iter++ )
    {
    /*
     *  update alpha
     */
    for(int t = 0; t < n_thread; t++ )
      dset( 2*k_Zf, ws[t], 0 );

    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < n_gene; j++ )
      {
      double *A = ws[omp_get_thread_num()];
      double *B = A + k_Zf;

      const uv *Y_j = Y + j*n_ac;
      uv *E_j = E + j*n_ac;
      for(int i = 0; i < n_ac; i++ )
        {
        if( w[i] == 0 ) continue;
        double W = w[i]*E_j[i].u/E_j[i].v;
        A[ Zf[i] ] += W*(Y_j[i].u - E_j[i].u);
        B[ Zf[i] ] += W*E_j[i].u;
        }
      }
    for(int t = 1; t < n_thread; t++ )
      daccum( 2*k_Zf, ws[0], ws[t] );

    // solve
    double *A = ws[0], *B = A + k_Zf;
    double mean = 0;
    for(int l = 0; l < k_Zf; l++ )
      {
      if(B[l] <= inv_tol ) A[l] = 0;
      else A[l] /= B[l];
      mean += A[l];
      }
    mean /= k_Zf;
    
    // constrain such that the mean is zero
    for(int l = 0; l < k_Zf; l++ )
      A[l] -= mean;

    // backtracking line-search
    for(int bt = 0; bt < bt_max; bt++ )
      {
      int iflag = 0;
      for(int l = 0; l < k_Zf; l++ )
        {
        B[l] = exp(A[l]);
        if(!isfinite(B[l])) iflag=1;
        }
      if(iflag) goto HALFSTEP;
      
      double dQ = 0;
      #pragma omp parallel for schedule(runtime) reduction(+:dQ)
      for(int j = 0; j < n_gene; j++ )
        {
        uv *E_j = E + j*n_ac;
        const uv *Y_j = Y + j*n_ac;
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i]==0) continue;
          double d = Y_j[i].u - E_j[i].u;
          double Dd2 = d*d;
          d = Y_j[i].u - E_j[i].u * B[Zf[i]];
          Dd2 -= d*d;
          dQ += w[i]*Dd2/E_j[i].v;
          }
        }
      if( dQ > bt_tol * n_ac * n_gene ) break;

      HALFSTEP:
      for(int l = 0; l < k_Zf; l++ )
        A[l] /= 2;
      nbt++;
      } 

    for(int l = 0; l < k_Zf; l++ )
      alpha[l] += A[l]; 

    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < n_gene; j++ )
      {
      uv *E_j = E + j*n_ac;
      for(int i = 0; i < n_ac; i++ )
        E_j[i].u *= B[Zf[i]];
      }

    /*
     * update beta and phi
     */
    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j < n_gene; j++ )
      {
      double *A = ws[ omp_get_thread_num() ];
      double *B = A + k_X1;

      const uv *Y_j = Y + j*n_ac;
      uv *E_j = E + j*n_ac;
      double *betaX1_j = beta + j*sz_beta; 

      // beta of X1
      dset( 2*k_X1, A, 0 );
      for(int i = 0; i < n_ac; i++ )
        {
        if(w[i]==0) continue;
        double W = w[i]*E_j[i].u/E_j[i].v;
        A[X1[i]] += W*(Y_j[i].u - E_j[i].u);
        B[X1[i]] += W*E_j[i].u;
        }
      for(int l = 0; l < k_X1; l++ )
        {
        if(B[l] <= inv_tol) A[l] = 0;
        else A[l] /= B[l];
        }

      for(int bt = 0; bt < bt_max; bt++ )
        {
        for(int l = 0; l < k_X1; l++ )
          while( !isfinite( B[l] = exp(A[l]) ) )
            A[l] /= 2;

        double dQ = 0;
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i]==0) continue;
          double d = Y_j[i].u - E_j[i].u;
          double Dd2 = d*d;
          d = Y_j[i].u - E_j[i].u * B[X1[i]];
          Dd2 -= d*d;
          dQ += w[i]*Dd2/E_j[i].v;
          }
        if( dQ > bt_tol * n_ac ) break;

        for(int l = 0; l < k_X1; l++ )
          A[l] /= 2;
        nbt++;
        }

      for(int l = 0; l < k_X1; l++ )
        betaX1_j[l] += A[l];
      for(int i = 0; i < n_ac; i++ )
        {
        if(w[i]==0 ) continue;
        E_j[i].u *= B[X1[i]];
        }

      // beta of Xf 
      double *betaXf_j = betaX1_j + k_X1;
      for(int p = 0; p < n_Xf; p++ )
        {
        const int *Xfp = Xf + p*n_ac;
        const int k_Xfp = k_Xf_[p];
        double *betaXf_jp = betaXf_j;

        dset( 2*k_Xfp, A, 0 );
        B = A + k_Xfp;
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i]==0) continue;
          double W = w[i]*E_j[i].u/E_j[i].v;
          A[Xfp[i]] += W*(Y_j[i].u - E_j[i].u);
          B[Xfp[i]] += W*E_j[i].u;
          }
        double mean = 0;
        for(int l = 0; l < k_Xfp; l++ )
          {
          if( B[l] == 0 ) A[l] = 0;
          else A[l] /= B[l];
          mean += A[l];
          }
        mean /= k_Xfp;
        for(int l = 0; l < k_Xfp; l++ )
          A[l] -= mean;

        for(int bt = 0; bt < bt_max; bt++ )
          {
          int iflag = 0;
          for(int l = 0; l < k_Xfp; l++ )
            {
            B[l] = exp(A[l]);
            if(!isfinite(B[l])) iflag =1;
            }
          if(iflag) goto HALFSTEP2;
          
          double dQ = 0;
          for(int i = 0; i < n_ac; i++ )
            {
            if(w[i]==0) continue;
            double d = Y_j[i].u - E_j[i].u;
            double Dd2 = d*d;
            d = Y_j[i].u - E_j[i].u * B[Xfp[i]];
            Dd2 -= d*d;
            dQ += w[i]*Dd2/E_j[i].v;
            }
          if(dQ > bt_tol * n_ac) break;

         HALFSTEP2:
          for(int l = 0; l < k_Xfp; l++ )
            A[l] /= 2;
          nbt++;
          }
        
        for(int l = 0; l < k_Xfp; l++ )
          betaXf_j[l] += A[l];
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i]==0) continue;
          E_j[i].u *= B[Xfp[i]];
          }
        
        betaXf_j += k_Xfp;
        }

      // beta of Xc
      double *betaXc_j = betaXf_j + k_Xf;
      for(int p = 0; p < n_Xc; p++ )
        {
        const double *Xcp = Xc + p * n_ac;
        double A = 0, B = 0;
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i]==0) continue;
          double W = w[i]*E_j[i].u/E_j[i].v * Xcp[i];
          A += W*(Y_j[i].u-E_j[i].u);
          B += W*E_j[i].u * Xcp[i];
          }
        if( B < inv_tol ) A = 0;
        else A /= B;
        while( !isfinite(B = exp(A)) )
          A /= 2;
        double dQ = 0;
        for(int bt = 0; bt < bt_max; bt++ )
          {
          for(int i = 0; i < n_ac; i++ )
            {
            if(w[i]==0) continue;
            double d = Y_j[i].u - E_j[i].u;
            double Dd2 = d*d;
            d = Y_j[i].u - E_j[i].u * B;
            Dd2 -= d*d;
            dQ += w[i]*Dd2/E_j[i].v;
            }
          if( dQ > bt_tol * n_ac ) break;
          while( !isfinite(B = exp(A))) A /= 2;
          nbt++;
          }
        betaXc_j[p] = A;
        for(int i = 0; i < n_ac; i++ )
          {
          if(w[i]==0) continue;
          E_j[i].u *= B;
          }
        }
      
      /*
       * update  phi
       */
      double *phi_j = phi + j*k_G1;
      B = A + k_G1;
      dset( 2*k_G1, A, 0 );

      for(int i = 0; i < n_ac; i++ )
        {
        if( w[i] == 0 ) continue;
        double W = w[i]*E_j[i].u * E_j[i].u / E_j[i].v;
        double dd = Y_j[i].u - E_j[i].u;
        dd *= dd;
        A[G1[i]] += W*dd /E_j[i].v;
        B[G1[i]] += W;
        }

      // backtracking doesn't seem to be necessary
      for(int l = 0; l < k_G1; l++ )
        if(B[l] > inv_tol && isfinite(A[l])) phi_j[l] *= A[l]/B[l];

      for(int i = 0; i  < n_ac; i++ )
        {
        if(w[i]==0) continue;
        E_j[i].v = Y_j[i].v + phi_j[G1[i]]*E_j[i].u*E_j[i].u;
        }
      }

    pll( &n_ac, &n_gene, w, Y, E, &Qd, &Qe );
    double Q = Qd + Qe;
    double progress = fabs(Q-Qold)/fabs(Qold);
    
    if(verbose == 1 )
      mess("*");
    else if(verbose > 1) 
      mess("[%4d]: %g %g\n", iter, Q/(n_ac*n_gene), progress);
    if( progress < iter_tol ) break;
    Qold = Q;
    }
  if(verbose==1) mess("\n");
  status_[0] = (iter < iter_max ? 0 : 1);
  status_[1] = iter;
  status_[2] = nbt;
  
  for(int t = 0; t < n_thread; t++ )
    free(ws[t]);
  return;
}
