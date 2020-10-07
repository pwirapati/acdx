#include <stdlib.h>
#include <string.h>

#include "etc.h"

spar *
spar_alloc(size_t szi, size_t szp)
{
  spar *S = malloc(sizeof(*S));
  S->i = malloc( sizeof(*S->i) * szi );
  S->x = malloc( sizeof(*S->x) * szi );
  S->szi = szi;
  S->pi = 0;
  S->p = malloc( sizeof(*S->p) * szp );
  S->szp = szp;
  S->p[0] = 0;
  S->pp = 0;
  return S;
}

void
spar_free(spar *S)
{
  free(S->i);
  free(S->x);
  free(S->p);
  free(S);
}

int
spar_addi( spar *S, int i, double x )
{
  if(S->pi >= S->szi)
    {
    void *t = realloc( S->i, sizeof(*S->i) * (S->szi * 2));
    if(t==NULL) return -2;
    S->i = t;

    t = realloc( S->x, sizeof(*S->x) * (S->szi * 2));
    if(t==NULL) return -2;
    S->x = t;
    
    S->szi *= 2;
    }
  
  S->i[S->pi] = i;
  S->x[S->pi] = x;
  S->pi++;
}

int
spar_addp( spar *S )
{
  if( S->pp >= S->szp )
    {
    void *t = realloc( S->p, sizeof(*S->p) * (S->szp*2));
    if(t==NULL) return -2;
    S->p = t;
    S->szp *= 2;
    }

  S->pp++;
  S->p[S->pp] = S->pi;
}

spar *
spar_transpose(spar *X, size_t n)
{
  spar *Xt = spar_alloc( X->pi, n + 1 );
  if( NULL == Xt ) return NULL;
  Xt->pi = X->pi;
  Xt->pp = n+1;

  for(size_t i = 0; i <= n; i++ )
    Xt->p[i] = 0;
  for(size_t k = 0; k < X->pi; k++ )
    Xt->p[ X->i[k] + 1 ]++;
  for(size_t i = 1; i <= n; i++ )
    Xt->p[i] += Xt->p[i-1];
  for(size_t j = 0; j < X->pp; j++ )
    for(size_t k = X->p[j]; k < X->p[j+1]; k++ )
      {
      size_t i = X->i[k];
      Xt->i[ Xt->p[i] ] = j;
      Xt->x[ Xt->p[i] ] = X->x[k];
      Xt->p[i]++;
      }
  for(size_t i = n-1; i > 0; i-- )
    Xt->p[i] = Xt->p[i-1];
  Xt->p[0] = 0;
  return Xt;
}

#define STRVEC_SZ_MIN (8)

svec *
svec_alloc(size_t szchr, size_t szvec)
{
  svec *S = malloc(sizeof(*S));

  if(szchr < STRVEC_SZ_MIN ) szchr = STRVEC_SZ_MIN;
  S->chr = malloc(sizeof(*(S->chr))*szchr);
  S->szchr = szchr;
  S->pchr = 0;

  if(szvec < STRVEC_SZ_MIN ) szvec = STRVEC_SZ_MIN;
  S->vec = malloc(sizeof(*(S->vec))*szvec);
  S->szvec = szvec;
  S->pvec = 0;
  S->vec[0] = 0;
  return S;
}

size_t
svec_append( svec *S, char *s )
{
  int len = strlen(s) + 1;
  while( S->pchr + len >= S->szchr )
    {
    char *chr = realloc(S->chr, sizeof(*(S->chr) ) * (S->szchr * 2));
    if(chr == NULL) return -2;
    S->szchr *= 2;
    S->chr = chr;
    }

  size_t i = S->pvec;
  strcpy( S->chr +  S->vec[i], s );
  S->pchr += len;
  
  S->pvec++;
  if( S->pvec >= S->szvec )
    {
    size_t *vec = realloc( S->vec, sizeof(*(S->vec)) * (S->szvec * 2));
    if(vec == NULL ) return -2;
    S->szvec *= 2;
    S->vec = vec;
    }
  S->vec[S->pvec] = S->pchr;

  return i;
}

void
svec_free(svec *S)
{
  free(S->chr);
  free(S->vec);
  free(S);
}
