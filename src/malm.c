#include "malm.h"
#include "nchol.h"

// pseudo-log-likelihood
double
PLL ( 
  const int n,
  const double *w,
  const uv *Y,
  const uv *E
  )
{
  double Ltry = 0;

  for(int i = 0; i < n; i++ )
    {
    if(w[i] == 0 ) continue;
    double d = Y[i].u - E[i].u;
    Ltry -= w[i]*d*d/E[i].v;
    if(isfinite(E[i].v)) Ltry -= w[i]*log(2*M_PI*E[i].v);
    }
return Ltry/2;
}

double
malm_init(
  const int n,
  const int p,
  const int q,
  const double *w,
  const uv *y,
  const double *Xt,
  const int *Gi,
  uv *E,
  const double *offset,
  double *beta,
  double *phi,
  double *S,
  double *H,
  int nnz,
  double d_adj
  )
{
  dset(p+p*(p+1)/2, S, 0);
  for(int i = 0; i < n; i++ )
    {
    if(w[i]==0) continue;
    double *H_b = H;
    const double *Xti = Xt + i*p;
    for(int b = 0; b < p; b++, H_b += b )
      {
      double W = w[i] * y[i].u * y[i].u / y[i].v;
      S[b] += W * Xti[b] * (log(y[i].u)-offset[i]);
      for(int a = 0; a <= b; a++ )
        H_b[a] += W * Xti[a] * Xti[b];
      }
    }
  nchol( p, H, NCHOL_AUTOTOL );
  nchol_fsub( p, H, S );
  nchol_bsub( p, H, S );
  dcp(p, S, beta );
  double Swrr[q];
  for(int a = 0; a < q; a++ ) Swrr[a]=0;
  for(int i = 0; i < n; i++ )
    {
    if(w[i]==0) continue;
    const double *Xti = Xt + i*p;
    double eta_i = 0;
    for(int b = 0; b < p; b++ )
      eta_i += Xti[b] * beta[b];
    E[i].u = exp(eta_i);
    double r = y[i].u / E[i].u - 1;
    Swrr[Gi[i]] += r*r;
    }
  for(int a = 0; a < q; a++ )
    phi[a] = Swrr[a]/nnz * d_adj;
  for(int i = 0; i < n; i++ )
    {
    if(w[i]==0) continue;
    E[i].v = y[i].v + phi[Gi[i]] * E[i].u * E[i].u;
    }

  return PLL(n,w,y,E);
}

void
update_phi( 
  int n,
  int q, 
  const double *w,
  const int *Gi,
  double *phi,
  const uv *y,
  uv *E,
  double d_adj
  )
{
  double Sd[q], Ss2[q];
  for(int a = 0; a < q; a++ )
    Sd[a] = Ss2[a] = 0;
  for(int i = 0; i < n; i++ )
    {
    if(w[i]==0)continue;
    double d = y[i].u - E[i].u;
    d = d*d * d_adj;
    double W = E[i].u/E[i].v;
    W *= w[i] * W;
    Sd[Gi[i]] += W * d;
    Ss2[Gi[i]] += W * E[i].v;
    }

  for(int a = 0; a < q; a++ )
    if( Ss2[a] > 0 ) phi[a] *= Sd[a]/Ss2[a];

  for(int i = 0; i < n; i++ )
    {
    if(w[i]==0) continue;
    E[i].v = y[i].v + phi[Gi[i]] * E[i].u * E[i].u;
    }

}

void
malm (
  const int *dims_,
  const uv *Y,
  const double *w,
  const double *Xt,
  const int *Gi,
  uv *Em,
  const double *offset,
  double *Beta,
  double *Phi,
  double *L_,
  const int *iopt_,
  const double *dopt_,
  int *n_iter
  )
{
  const int n = dims_[0], m = dims_[1], p = dims_[2], q = dims_[3];

  const int verbose = iopt_[0];
  const int iter_max = iopt_[1];
  const int bt_iter_max = iopt_[2];

  const double epsilon = dopt_[0];
  const double bt_tol = dopt_[1];
  const double bt_step = dopt_[2];

  double *tmp = (double*)malloc( sizeof(double) * (p+p*(p+1)/2));

  double *S = tmp;
  double *H = S + p;
  int nnz = 0;
  for(int i = 0; i < n; i++ )
    if(w[i] > 0) nnz++;
  double d_adj = (double)nnz/(nnz-p); // correction for estimating the mean params

  for(int j = 0; j < m; j++ )
    {
    const uv *y = Y + j*n;
    uv *E = Em + j*n;
    double *beta = Beta + j*p;
    double *phi = Phi + j*q;

    double L = malm_init(n,p,q,w,y,Xt,Gi,E,offset,beta,phi,S,H,nnz,d_adj);
    double Ltry;

    int betaconv = 0, phiconv = 0;
    for(int iter = 1; iter <= iter_max; iter++ )
      {
      dset( p+p*(p+1)/2, S, 0);
      for(int i = 0; i < n; i++ )
        {
        if(w[i]==0) continue;
        double *H_b = H;
        const double *Xti = Xt + i*p;
        for(int b = 0; b < p; b++, H_b += b )
          {
          double W = w[i] * E[i].u * E[i].u / E[i].v;
          S[b] += W * Xti[b] * (y[i].u/E[i].u - 1);
          for(int a = 0; a <= b; a++ )
            H_b[a] += W * Xti[a] * Xti[b];
          }
        }
      nchol( p, H, NCHOL_AUTOTOL );
      nchol_fsub( p, H, S );
      nchol_bsub( p, H, S );

      int bt = 0;

    BACKTRACK:
      for(int i = 0; i < n; i++ )
        {
        const double *Xti = Xt + i*p;
        double eta_i = offset[i];
        for(int a = 0; a < p; a++ )
          eta_i += Xti[a]*(beta[a] + S[a]);
        E[i].u = exp(eta_i);
        E[i].v = y[i].v + phi[Gi[i]] * E[i].u * E[i].u;
        }
      
      double phi_t[q];
      dcp( q, phi, phi_t);
      update_phi(n,q,w,Gi,phi_t,y,E,d_adj);

      Ltry = PLL(n,w,y,E);

      if( bt < bt_iter_max && (Ltry - L)/(nnz-p-q) < - bt_tol ) // backtrack
        {
        for(int a = 0; a < p; a++ )
          S[a] *= bt_step;
        bt++;
        goto BACKTRACK;
        }
      for(int a = 0; a < p; a++ )
        beta[a] += S[a];
      dcp( q, phi_t, phi );
      
      if(verbose == 2) mess("[%3d] %g (%+g)\n",iter,L, Ltry-L);
      if( fabs(Ltry - L)/(nnz - p - q) < epsilon )
        {
        if(n_iter) n_iter[j] = iter;
        break;
        }
      L = Ltry;
      } // iter
    
    if(L_) L_[j] = L;
    
    // flag dropped terms
    for(int a = 0; a < p; a++ )
      if( H[a*(a+1)/2-1] == 0) beta[a] = NAN;

    } // gene

  free(tmp);
  return;   
}

