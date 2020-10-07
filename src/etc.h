#ifndef _ETC_H_
#define _ETC_H

typedef struct
  {
  int *i;
  double *x;
  size_t szi;
  size_t pi;

  int *p;
  size_t szp;
  size_t pp;
  }
spar;

spar* spar_alloc(size_t szi, size_t szp);
void spar_free(spar *S);
int spar_addi( spar *S, int i, double x );
int spar_addp( spar *S );
spar* spar_transpose(spar *X, size_t n);

typedef struct 
  {
  char *chr;
  size_t szchr;
  size_t pchr;

  size_t *vec;
  size_t szvec;
  size_t pvec;
  }
svec;

svec* svec_alloc(size_t szchr, size_t szvec);
size_t svec_append( svec *S, char *s );
void svec_free(svec *S);

#endif // _ETC_H_
