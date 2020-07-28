#include <stdio.h>
#include "nglm.h"

void
R_nglm_gamma_ss(
  int *dims,
  double *w, double *Y, int *Xp, int* Zq,
  double *alpha, double *v_alpha, double* beta, double *v_beta,
  double *phi,  double *E,
  double *shift,
  double *tol, int *itermax
  )
{
  int ni = dims[0], nj = dims[1], np = dims[2], nq = dims[3];
  nglm_gamma_ss(
    ni, nj, np, nq,
    w, (uv*)Y, Xp, Zq,
    alpha, v_alpha, beta, v_beta, phi, (uv*)E,
    shift[0], 
    tol[0],itermax[0]
    );
}
void
R_nglm_gamma_sd(
  int *dims,
  double *w, double *Y, double *Xt, int* Zq,
  double *alpha, double *v_alpha, double* beta, double *v_beta,
  double *phi,  double *E,
  double *shift,
  double *tol, int *itermax
  )
{
  int ni = dims[0], nj = dims[1], np = dims[2], nq = dims[3];
  nglm_gamma_sd(
    ni, nj, np, nq,
    w, (uv*)Y, Xt, Zq,
    alpha, v_alpha, beta, v_beta, phi, (uv*)E,
    shift[0], 
    tol[0],itermax[0]
    );
}
