/*
 * tab2spar.c - read tab-delimited data table from R connection into
 *              a sparse `dgCMatrix` format
 *
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Connections.h>

#include "etc.h"

typedef struct
  {
  Rconnection conn;
  char *buf;
  size_t szbuf;
  size_t current;
  size_t last;
  }
wcon;

wcon*
wcon_alloc( Rconnection conn, size_t szbuf )
{
  wcon *F = (wcon*) malloc(sizeof(wcon));
  F->conn = conn;
  F->szbuf = szbuf;
  F->buf = (char*)malloc(sizeof(char) * szbuf );
  F->current = 0;
  F->last = 0;
  return F;
}

void
wcon_free( wcon *F)
{
  free(F->buf);
  free(F);
}

#define MALLOC_ERROR (-2)
#define EOFILE (-1)
#define EOCOL (0)
#define EOROW (1)

int
get_str( wcon *F, char **s_)
{
  int status;
  size_t i = F->current; 
  while(1)
    {
    if( i >= F->last )  // reach end of buffer, need to get more
      {
      i -= F->current; // i is the length of the current (partial) word
      if( i > F->szbuf/2 )
        {
        size_t newszbuf = F->szbuf * 2;
        char *tmp = (char*) realloc(F->buf, newszbuf );
        if(tmp == NULL) return MALLOC_ERROR;
        F->buf = tmp;
        F->szbuf = newszbuf;
        }
      memcpy(F->buf, F->buf + F->current, sizeof(char)*i);

      size_t n_read = R_ReadConnection(F->conn, F->buf + i, F->szbuf - i - 1);
      if(n_read == 0 ) 
        { status = EOFILE; break; }
      F->last = i + n_read;
      F->current = 0;
      }
    if( F->buf[i] == '\t' )
      { status = EOCOL; break; }
    else if( F->buf[i] == '\r' || F->buf[i] == '\n' )
      { status = EOROW; break; }
    else 
      i++;
    }
  F->buf[i] = '\0';

  *s_ = F->buf + F->current;
  F->current = i+1;
  return status;
} 

SEXP tab2spar(SEXP conn, SEXP transpose )
{
 
  wcon *F = wcon_alloc( R_GetConnection(conn), 64 );
 
  svec *H = svec_alloc(32*10000, 10000);

  char *s;
  int ncol = 0, r;
  get_str(F, &s); // discard header of first column
  do
    {
    r = get_str(F, &s);
    if( r < 0 ) break;
    if( -2 == svec_append( H, s ))
      fprintf(stderr,"realloc error");
    ncol++;
    }
  while( r == EOCOL );

  int nz = 0, nij = 0;
  svec *G = svec_alloc(32*10000,10000 );
  spar *X = spar_alloc(10000,100);

  int nrow = 0;
  while(1)
    {
    r = get_str(F, &s );
    if( r == EOFILE ) break;
    svec_append( G, s );

    int j = 0;
    do
      {
      r = get_str(F, &s);
      if( r < 0 ) break;
      double v = strtod( s, NULL );
      if(v != 0 )
        spar_addi( X, j, v);
      j++;
      }
    while( r == EOCOL );
    nrow++;
    spar_addp( X );
    if(nrow % 100 == 0 )
      fprintf(stderr,"\r%d rows read", nrow );
    }
  fprintf(stderr,"\r%d rows read\n",nrow);

  wcon_free(F);

  SEXP out = R_NilValue;
  if( X->pi > (size_t)INT_MAX )
    {
    fprintf(stderr,"sparse data has larger than %d non zeroes\n",INT_MAX);
    goto CLEANUP;
    }

  if( !LOGICAL(transpose)[0] )
    {
    fprintf(stderr,"transposing...\n");
    spar *Xt = spar_transpose(X, ncol);
    spar_free(X);
    X = Xt;

    svec *T = H; H = G; G = T;

    int t = nrow; nrow = ncol; ncol = t;
    }

  fprintf(stderr,"packaging to R objects...\n");
  SEXP colnames = PROTECT( allocVector(STRSXP,ncol) );
  for(int i = 0; i < ncol; i++ )
    SET_STRING_ELT(colnames,i,mkChar( H->chr + H->vec[i] ));

  SEXP rownames = PROTECT( allocVector(STRSXP,nrow) );
  for(int i = 0; i < nrow; i++ )
    SET_STRING_ELT(rownames,i,mkChar( G->chr + G->vec[i] ));

  out = PROTECT( allocVector(VECSXP,5));
  SET_VECTOR_ELT(out,0,colnames);
  SET_VECTOR_ELT(out,1,rownames);

  SEXP icol = PROTECT( allocVector( INTSXP, X->pi ));
  memcpy( INTEGER(icol), X->i, sizeof(*X->i) * X->pi);
  SET_VECTOR_ELT(out,2,icol);

  SEXP x = PROTECT( allocVector( REALSXP, X->pi ));
  memcpy( REAL(x), X->x, sizeof(*X->x) * X->pi);
  SET_VECTOR_ELT(out,3,x);
  
  SEXP prow = PROTECT( allocVector( INTSXP, X->pp + 1 ));
  memcpy( INTEGER(prow), X->p, sizeof(*X->p) * X->pp);
  INTEGER(prow)[X->pp] = X->pi;
  SET_VECTOR_ELT(out,4,prow);

  UNPROTECT(6);

CLEANUP:

  svec_free(H);
  svec_free(G);
  spar_free(X);

  fprintf(stderr,"Done.\n");

  return out;
}
