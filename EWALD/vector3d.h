#define DIM 3

#ifdef VECTORCACHE
typedef real vector[4];
#else
typedef real vector[DIM];
#endif

/* NOTE: 
  typedef real vector[DIM+1]; 
gives better fit of cache line and by about 1% faster performance of cook
However, polarizability does not work in this way ...
*/

#define VO(A,O) { A[0] O; A[1] O; A[2] O; }
#define VV(A,B) { A[0] B[0]; A[1] B[1]; A[2] B[2]; }
#define VVO(A,B,O) { A[0] B[0] O; A[1] B[1] O; A[2] B[2] O; }
#define VVV(A,B,C) { A[0] B[0] C[0]; A[1] B[1] C[1]; A[2] B[2] C[2]; }
#define VVVO(A,B,C,O) { A[0] B[0] C[0] O; A[1] B[1] C[1] O; A[2] B[2] C[2] O; }
#define VVVV(A,B,C,D) \
  { A[0] B[0] C[0] D[0]; A[1] B[1] C[1] D[1]; A[2] B[2] C[2] D[2]; }
#define VVVVV(A,B,C,D,E)                                                 \
  { A[0] B[0] C[0] D[0] E[0]; A[1] B[1] C[1] D[1] E[1]; A[2] B[2] C[2] D[2] E[2]; }
#define SQR(A) (A[0]*A[0]+A[1]*A[1]+A[2]*A[2])
#define SQRD(A,B) (Sqr(A[0]-B[0])+Sqr(A[1]-B[1])+Sqr(A[2]-B[2]))
#define SCAL(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define SUM(A) (A[0]+A[1]+A[2])
#define PROD(A) (A[0]*A[1]*A[2])
#define VECT(A,B,C) { \
  A[0]=B[1]*C[2]-B[2]*C[1]; \
  A[1]=B[2]*C[0]-B[0]*C[2]; \
  A[2]=B[0]*C[1]-B[1]*C[0]; }

#define VARG(A) (double)A[0],(double)A[1],(double)A[2]
#define AVARG(A) A,A+1,A+2
