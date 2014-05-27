#ifdef __cplusplus
#define TYPEOF(X) (typeof(X))
extern "C" {
#else
#define TYPEOF(X) /* no cast needed */
#endif

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "loop.h"

#define PI_DOUBLE 3.14159265358979323846
#ifndef PI
#define PI PI_DOUBLE
#endif

#ifndef Max
/* A=max(A,B) */
#define Max(A,B) { if ((B)>(A)) (A)=(B); }
/* A=min(A,B) */
#define Min(A,B) { if ((B)<(A)) (A)=(B); }
#endif

#ifndef max
/* see below for fmax and fmin */
#define max(A,B) ((A)<(B)?(B):(A))
#define min(A,B) ((A)<(B)?(A):(B))
#endif

double fmin(double x,double y) { if (x<y) return x; else return y; }
double fmax(double x,double y) { if (x>y) return x; else return y; }

#ifndef copy
#ifdef __TRANSPUTER__
#define copy block_move
void block_move(void *, const void *, int);
#else
#define copy memcpy
#endif
#endif

#define Sqr(X) ((X)*(X))
#define Cub(X) ((X)*(X)*(X))
#define Pow4(X) Sqr(Sqr(X))
#define Pow5(X) ((X)*Sqr(Sqr(X)))
#define Pow6(X) Cub(Sqr(X))

// cbrt(): C99, POSIX.1-2001
#define cubrt(X) cbrt(X)
//#define cubrt(X) pow((X),(1e0/3))

#define Sign(X) ((X)<0?-1:(X)>0?1:0)

#define put(_X) printf("%11s=%-13.6g\n",#_X,(double)(_X));
#define put_(_X) printf("%11s=%-13.6g ",#_X,(double)(_X));
#define put2(_X,_Y) { put_(_X) put(_Y) }
#define put3(_X,_Y,_Z) { put_(_X) put_(_Y) put(_Z) }
#define _n putc('\n',stdout);
#define putv(_X) prt("%11s=(%13.6g %13.6g %13.6g )", \
  #_X,(double)_X[0],(double)_X[1],(double)_X[2]);

#define alloc(X,S) do { if (!(X=TYPEOF(X)malloc(S))) { fprintf(stderr,"no heap for " #X "\n"); exit(1); } } while (0)
#define alloczero(X,S) do { if (!(X=TYPEOF(X)malloc(S))) { fprintf(stderr,"no heap for " #X "\n"); exit(1); } memset(X,0,S); } while (0)

#define Error(X) do { fprintf(stderr,"ERROR %s\n",X); exit(1); } while(0)

/*.....#ifdef __sun__*/
/*.....#include "missing.h"*/
/*.....#endif*/

#include <stdarg.h>
#include "mystring.c"
#define alloc2Darray(X,N,M) do { int _I;   \
  alloc((X),(N)*sizeof((X)[0]));           \
  alloc((X)[0],(N)*(M)*sizeof((X)[0][0])); \
  for (_I=1; _I<(N); _I++) (X)[_I]=(X)[0]+_I*(M); } while (0);
#define alloc2Darrayzero(X,N,M) do { int _I;   \
  alloc((X),(N)*sizeof((X)[0]));           \
  alloczero((X)[0],(N)*(M)*sizeof((X)[0][0])); \
  for (_I=1; _I<(N); _I++) (X)[_I]=(X)[0]+_I*(M); } while (0);
#define free2Darray(X) do { free(X[0]); free(X); } while (0);

/* added 3/2004 */
#define arrayzero(X,Y) memset((X),0,(Y)*sizeof((X)[0]))
/* WARNING: this is not ANSI C portable, although I do not know about
   architecture where this trick would fail */
/* added 4/2008 (was already in alloc.h) */
#define copyarray(X,Y,Z) copy((X),(Y),(Z)*sizeof(Y[0]))

#define sdsalloc(X,Y) { alloc((X),Y); *(int*)(X)=Y; }

  /* MISSING \n !!! */
#define prt printf
#define ERROR(_X) { prt _X; prt("\n"); exit(0); }
#define WARNING(_X) { prt _X; prt("\n"); }

#ifdef __cplusplus
}
#endif

