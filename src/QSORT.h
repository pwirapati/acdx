/*
QSORT.h - macro for quicksort

@2010 Pratyaksha J. Wirapati

SYNOPSIS
--------
This is a macro to quicksort an array using a
comparison function defined by an arbitrary C expression passed
as an argument to the macro.

USAGE
-----

QSORT( itype, n, otype, o, CMP );

itype   type of array indices (integer types: [unsigned] char,short,int,long )
n       is the length of the sorted array
otype   type of array element to be sorted
o       array to be sorted, o[0] to o[n-1] will be affected 
CMP     an expression that contains two identifiers, _u and _v,
        which evaluates to one (or any non-zero) if _u "less than" _v.
        "less_than" means _u will appear to left of _v in the sorted list.
        The type of _u and _v is 'itype'.


Typically, indirect sorting is used for arbitrary data indexed by o.
It is possible to use this macro to directly sort o by using
  (_u < _v)
as the macro.

The index array 'o' does not have to contain 0..n-1. It may contain
already permuted indices, subset of full range of indices for the
sorted data, or even duplicated indices.

EXAMPLES
--------

1. Direct sorting

  QSORT( int, n, int, x, _u < _v )
  QSORT( int, n, double, x, _u < _v )
  QSORT( int, n, char*, x, strcmp(_u,_v) > 0 )   // descending order

2. Indirect sort of a double array in decreasing order
   (but keep the original ordering in case of ties):

  double *x, int *o;
  ...
  QSORT( int, n, int, o, 
    ( x[_u] > x[_v] ? 1 : (x[_u]==x[_v] ? (_u < _v) : 0) : 0 )
    );

  Note: quicksort can only be used to produce stable sorting
  when indirect.

3. Indirect sorting of a more complex data

  double *x, char **y, int *o;
  ...
  QSORT( int, n, int, o,
      ( c = strcmp( y[_u], y[_v] ), 
        (c < 0 ? 1 : ( c == 0 ? (x[_u] < x[_v] ? 1 : 0) : 0 ) )
      )
    );
        

TIPS
----

1. If the sorting with the same CMP, particularly complex ones,
   done repeatedly, it's better to define them separately as a function:

    void mysort(int n, int *o, complextype *data)
    {
      QSORT( int, n, int, o, CMP_USING_data )
    }

2. Inline function can be used to make CMP more readable.

   inline int myCMP(u,v,data)
    {
    return u_less_than_v
    }

   ...

   QSORT(n,o, myCMP(_u,_v,data) )

   However, if the compiler doesn't really inline the function, 
   repeated callbacks from the inner loop will be inefficient.
   Nevertheless it is still better than using qsort() from <stdlib.h>
   because we can pass the data in a thread-safe way.

CAUTIONS
--------
As usual with macros, n and o should not be expression with side effect
(such as those with post/pre-increment operator). 

The identifiers
  _itype, _n, _otype, _o, _CMP,
  _s, _k, _i, _j, _a, _b, _c, _t, _u, and _v
are locally defined inside the macro.
Except for _u and _v, they should never be used inside CMP.
The definition of _u and _v overshadows the same identifiers (if any)
in the calling context.

*/

#ifndef _QSORT_H_
#define _QSORT_H_

#ifndef QSORT_BLOCK
#define QSORT_BLOCK 32    // partition size handled by insertion sort
#endif

#define QSORT( _itype, _n, _otype, _o, _CMP )                                \
{                                                                            \
  int _s[ 8*sizeof(_otype) ];                                                \
  int _k = 0;  /* stack pointer */                                           \
                                                                             \
  _itype _i, _j, _a, _b, _c;                                                 \
  _otype _t;                                                                 \
  _otype _u;                                                                 \
  _otype _v;                                                                 \
  _a = 0; _b = _n-1;                                                         \
  for(;;)                                                                    \
    {                                                                        \
    while( _b - _a >= QSORT_BLOCK )                                          \
      {                                                                      \
      /* median of three */                                                  \
      _c = ( _a + _b )/2;                                                    \
      if( _u = _o[_b], _v = _o[_a], _CMP )                                   \
        { _t = _o[_a]; _o[_a] = _o[_b]; _o[_b] = _t; }                       \
      if( _u = _o[_c], _v = _o[_b], _CMP )                                   \
        {                                                                    \
        _t = _o[_c]; _o[_c] = _o[_b];                                        \
        if( _u = _t, _v = _o[_a], _CMP )                                     \
          { _o[_b] = _o[_a]; _o[_a] = _t; }                                  \
        else                                                                 \
          _o[_b] = _t;                                                       \
        }                                                                    \
                                                                             \
      /* main exchange loop */                                               \
      _i = _a-1; _j = _b;                                                    \
      for(;;)                                                                \
        {                                                                    \
        while( _u = _o[++_i], _v = _o[_b], _CMP ) ;                          \
        while( _u = _o[_b], _v = _o[--_j], _CMP ) ;                          \
                                                                             \
        if( _i >= _j ) break;                                                \
        _t = _o[_i]; _o[_i] = _o[_j]; _o[_j] = _t;                           \
        }                                                                    \
      _t = _o[_i]; _o[_i] = _o[_b]; _o[_b] = _t;                             \
                                                                             \
      if( _i - _a > _b - _i )                                                \
        { _s[_k++] = _a; _s[_k++] = _i - 1; _a = _i + 1; }                   \
      else                                                                   \
        { _s[_k++] = _i + 1; _s[_k++] = _b; _b = _i - 1; }                   \
      }                                                                      \
    if( _k == 0 ) break;                                                     \
    _b = _s[--_k]; _a = _s[--_k];                                            \
    }                                                                        \
                                                                             \
  /* use insertion sort for the remaining small partitions */                \
  for( _i = 1; _i < _n; _i++ )                                               \
    {                                                                        \
    _t = _o[_i];                                                             \
    for( _j = _i; _j > 0 && (_u = _t, _v = _o[_j-1], _CMP); _j-- )           \
      _o[_j] = _o[_j - 1];                                                   \
    _o[_j] = _t;                                                             \
    }                                                                        \
}                                                                            \


#endif // _QSORT_H_ 
