#ifndef _NGLM_H_
#define _NGLM_H_

/*
 * nglm.h
 */


typedef struct {
  double u;
  double v;
  }
uv;

extern
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
;

extern
int
nglm_gamma_sd (
  const int ni,
  const int nj,
  const int nk,
  const int nl,

  const double *w,
  const uv *Y,
  const double *Xt,
  const int *Zl,

  double *alpha,          // 1*nl
  double *v_alpha,        // 1*1*nl
  double *beta,           // nk*nj
  double *v_beta,         // nk*nk*nj
  double *phi,            // 1*nj
  uv *E,                  // ni*nj
  
  const double shift,     // response = log( y + shift )
 
  const double tol,
  const int itermax
  )
;

#endif // _NGLM_H_
