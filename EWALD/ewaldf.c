/* cc -O2 -o ewaldf ewaldf.c -lm

This is a serial nonpolar simply-interfaced version of MACSIMUS/sim/ewald.c

This code uses ideas of J. Perram (making use of symmetry).
See also macsimus/c/ewald.c (fool-proof slow standalone version for debugging).

  pass -1: silent initialization
        0: initialization & protocol printed
        1: sums Q(k)
        2: energy (return value) and forces calculated

Calling scheme:
^^^^^^^^^^^^^^^
vector boxL; // simulation box
vector *forces,*configuration; // [el.No]
real *charges; // [el.No]
Ewald_k_space(0,<any>,boxL,<any>,<any>);  // initialization - may be repeated, see el
initerfc(grid, minr,maxr, el.cutoff, el.alpha, 1, 3); // if splines
repeat {
  forces=0;
  Ewald_k_space(1,<any>,boxL,charges,configuration);
  Ek=Ewald_k_space(2,forces,boxL,charges,configuration);
  Er=Ewald_r_space[_splines](forces,boxL,charges,configuration);
}

Formulas:
^^^^^^^^^
U =
    SUM_j<l qj ql Erfc(alpha rjl)/rjl [ r-space term: see module interpol.c ]
  + SUM_k!=0 exp(-pi^2/alpha^2*kk) / (2 V pi kk) |Q(k)|^2
  + 2 pi / (2 epsinf+1) *M^2/V
  - alpha/sqrt(pi)*SUM_j qj^2
(is divided by 4 pi eps0 in SI)
where
  kk denotes (kx/Lx)^2+(ky/Ly)^2+(ky/Ly)^2
  Q(k) = SUM_j qj exp[2 pi i (kx/Lx.xj+ky/Ly.yj+kz/Lz.zj)]
  M = SUM_j qj rj
  V = Lx Ly Lz

More:
^^^^^
For ready-to-code formulas, see Aguado, Madden: JCP 119, 7471 (2003).
Differences:
* h=2*pi*k/L (vector, by components)
* r-space parameter is called kappa (here: alpha, while kappa means k-space)
* constant term -alpha/sqrt(pi)*SUM_j qj^2 is missing in U

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#define PRESSURETENSOR 0
//#define PRESSURETENSOR (PT_VIR|PT_OFF)

typedef double real;
#include "include.h"
#include "vector3d.h"
#include "ewaldf.h"
struct el_s el;

/* global : */
Q_t *Q;

static double estimr,estimk;

/*** complex arithmetic ***/
#define CxE(A) { A.re=(c=A.re)*E->re-A.im*E->im; A.im=c*E->im+A.im*E->re; }
#define CxiE(A) { A.re=(c=A.re)*E->re+A.im*E->im; A.im=A.im*E->re-c*E->im; }
#define Cadd(A,B) { A.re+=B.re; A.im+=B.im; }
#define Csqr(A) (A.re*A.re+A.im*A.im)

/*** serial: loops over charges ***/
#define iqloop for (iqt=0; iqt<Nq; iqt++)
#define iqfloop for (iqt=0; iqt<Nq; iqt++)
#define iqeloop(E0) for (iqt=0; E=E0+iqt,iqt<Nq; iqt++)

typedef struct {
  complex
    x,   /* q*exp[2*PI*i*x*kx/Lx] */
    xy,  /* q*exp[2*PI*i*(x*kx/Lx+y*ky/Ly)] */
    xY,  /* q*exp[2*PI*i*(x*kx/Lx-y*ky/Ly)] */
    xyz, /* q*exp[2*PI*i*(x*kx/Lx+y*ky/Ly+z*kz/Lz)] */
    xyZ, /* q*exp[2*PI*i*(x*kx/Lx+y*ky/Ly-z*kz/Lz)] */
    xYz, /* q*exp[2*PI*i*(x*kx/Lx-y*ky/Ly+z*kz/Lz)] */
    xYZ; /* q*exp[2*PI*i*(x*kx/Lx-y*ky/Ly-z*kz/Lz)] */
  complex dummy; /* cache-line optimization */
} qtab_t;

/*** some ugly globals with sums of charges etc. ***/
struct charges_s charges;

#define QTIMES /* dummy here */

/* moved here because of PARALLEL */
static int Nq;               /* # of all charges */
static int kxmax;
static complex *ex,*ey,*ez;  /* array[Nq] : exp(2*PI*i*r/L) */
static qtab_t *qtab;

static real *ekk;    /* array[Nkk] : table of ...exp(...k*k) factors */
static struct ktab_s {
  int kymax;
  int *kzmax; /*[kymax(incl.)]*/
} *ktab /*[kxmax(incl.)]*/ ;
static int ktabsize; /* =kxmax+1 */
static vector *F;    /* array[Nq] */
static real machineprec;
static real KX,PIaq;

static int almostint(real x) /************************************ almostint */
/* returns true if x is very close to an integer */
{
  x-=(int)x;
  return fabs(0.5-fabs(x-0.5))<machineprec;
}

static int iscube(vector boxL) /************************************* iscube */
/* returns 1 if the box is a cube (with machine precision) */
{
  int i;

  loopto (i,1,2)
    if (fabs(boxL[0]-boxL[i])/(boxL[0]+boxL[i])>1e-7) return 0;

  return 1;
}

/*                                                     ****** Ewald_k_space */
double Ewald_k_space(int pass,vector *frp,vector boxL,vector *rp,real *q)
{
  static vector oldL={-3e33,-3e33,-3e33};
  static vector firstL;

  static int Nkk= -1; /* # of k-vectors irresp. of signs of components */
                      /* -1 forces initialization */

  real *ptr;
  vector ek;
  complex *Qk;
  real e,c,kk,aux,kxq,kyq;
  vector PIL;
  static real Kk; /* K refers to z-coordinate integer vector */
  real dK;
  double Energy=0,alsPI=el.alpha/sqrt(PI);
  int i,kx,ky,kz,iQ,ikk,kymax,kzmax;
  long size,qtabsize;
  real FA,FB,FC,FD;
  complex *E,A,B,C,D;
  int jkk,iqt;

  int iq,n,ns;
  vector *r,*fr;
  int nagain=0;
  static int oldkappa=-99;

  if (el.kappa<=0) return 0;

  VV(PIL,=2*PI/boxL)

  if (pass<=0) {
    /*
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %                        pass 0 - initialization                       %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    */

    if (Nkk==-1 || fabs(el.kappa-oldkappa)>1e-11) {
      oldkappa=el.kappa;
      if (machineprec==0) {
        if (sizeof(machineprec)<8) machineprec=1e-6;
        else machineprec=1e-9; }
      dK=machineprec*2;

     AGAIN: if (nagain++>25) ERROR(("Ewald: cannot solve rounding problems"))

      if (pass>=0) prt("initializing Ewald: kappa=K/L=%g\n",el.kappa);
      KX=boxL[0]*el.kappa;
      Kk=Sqr(el.kappa);

      /*** if not the 1st initialization then free the previous arrays ***/
      if (Nkk != -1) {
        if (ekk) free(ekk);
        if (Q) free(Q);
        if (ktab) {
          loop (n,0,ktabsize) free(ktab[n].kzmax);
          free(ktab); } }

      /*** computing Nq and (sums of) charges ***/
      Nq=0;
      charges.sum=charges.sumq=charges.max=0;

      loop (n,0,el.No) {
        if ((e=q[n])!=0) {
          Max(charges.max,fabs(e))
          charges.sum+=e; charges.sumq+=e*e;
          Nq++; }
      } /* n */

      /*** computing diagonal correction - REMOVED */

      /*** computing sizes of Q and ekk;
           readjusting K for rounding problems:
             because of optimization, an expression like
             L[2]*sqrt(kyq) may give once e.g. 2+1e-15 and in other part
             of the code 2-1e-15; then, (int) of such a number is not
             well defined: K is a bit increased if this might happen;
           (computing ekk moved to pass==1 with check of L change) ***/
      iQ=ikk=0;
      kxmax=(int)KX;
      loopto (kx,0,kxmax) {
        kxq=Kk-Sqr(kx/boxL[0]);
        aux=boxL[1]*sqrt(kxq);
        if (pass==0 && almostint(aux)) {
          el.kappa+=el.kappa*dK; dK*=1.4;
          prt("Ewald: possible rounding problem (x) => kappa increased to %.12g\n",el.kappa);
          goto AGAIN; }
        kymax=(int)aux;
        loopto (ky,0,kymax) {
          kyq=kxq-Sqr(ky/boxL[1]);
          aux=boxL[2]*sqrt(kyq);
          if (pass==0 && almostint(aux)) {
            el.kappa+=el.kappa*dK; dK*=1.4;
            prt("Ewald: possible rounding problem (y) => kappa increased to %.12g\n",el.kappa);
            goto AGAIN; }
          kzmax=(int)aux;
          loopto (kz,!(kx+ky),kzmax) {
            iQ += 1 << ((kx>0)+(ky>0)+(kz>0)-1); /* = 2^(..) */
            ikk++; } } } /* kz,ky,kx */

      /*** creating ktab (repeats the previous code to some extent) ***/
      ktabsize=kxmax+1;
      allocarray(ktab,ktabsize);
      loopto (kx,0,kxmax) {
        kxq=Kk-Sqr(kx/boxL[0]);
        kymax=(int)(boxL[1]*sqrt(kxq));
        ktab[kx].kymax=kymax;
        allocarray(ktab[kx].kzmax,kymax+1);
        loopto (ky,0,kymax) {
          kyq=kxq-Sqr(ky/boxL[1]);
          ktab[kx].kzmax[ky]=(int)(boxL[2]*sqrt(kyq)); } }

      VO(oldL,=-3e33) /* to pretend the old L is different and calculate ekk later */

      if (!pass) { /*** initialization protocol printed ***/
#include "ewaldprt.c"
      } /* !pass */

      sdsalloc(Q,sizeof(Q_t)-sizeof(complex)+ iQ*sizeof(complex))
      Q->N=iQ;
      VO(Q->M,=0) /* to have Q->M defined even if no charges
                     ... though it is stupid to allocate Q at all */
      Nkk=ikk;
      /* during pre-inititalization, Nkk may be zero, so prevent alloc error */
      alloc(ekk,(Nkk+!Nkk)*sizeof(real)); }

    return (double)Q->N/el.No;
  } /* initialization (pass<=0) */

  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             passes 1 and 2                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */

  /*** computing ekk (moved here from initialization to enable L change) ***/
  /* selected k-vectors, table sizes, alpha, ... are unchanged! */
//  if (L[0]!=oldL[0] || L[1]!=oldL[1] || L[2]!=oldL[2]) {
  if (memcmp(boxL,oldL,sizeof(boxL))) {
    if (firstL[0]!=0)
      if (fabs(log(boxL[0]/firstL[0]))>el.diff
       || fabs(log(boxL[1]/firstL[1]))>el.diff
       || fabs(log(boxL[2]/firstL[2]))>el.diff) {
        WARNING(("variable L with Ewald:\n*** box has changed by >%g%% (log-scale) and k-vectors remain unchanged\n*** L=%g %g %g",el.diff*100,VARG(boxL)))
        el.diff*=2; }

    if (oldL[0]>-1e33 && firstL[0]==0) {
      prt("WARNING: variable box L with Ewald:\n(only small changes in L are allowed without restarting the simulation)\n");
      VV(firstL,=oldL) }

    PIaq=-Sqr(PI/el.alpha);
    iQ=ikk=0;
    if (kxmax+1!=ktabsize)
      ERROR(("ktab size error: requested=%d  allocated with=%d  KX=%g",kxmax+1,ktabsize,KX))
    loopto (kx,0,kxmax) {
      kymax=ktab[kx].kymax;
      loopto (ky,0,kymax) {
        aux=Sqr(ky/boxL[1])+Sqr(kx/boxL[0]);
        kzmax=ktab[kx].kzmax[ky];
        loopto (kz,!(kx+ky),kzmax) {
          kk=Sqr(kz/boxL[2])+aux;
          ekk[ikk] = 4/(kk*(boxL[0]*boxL[1]*boxL[2]))*exp(PIaq*kk);
          /* NB: factor 4 above comes from: 2 from differentiation of Q^2
                                            2 from counting k and -k at once */
          iQ += 1 << ((kx>0)+(ky>0)+(kz>0)-1); /* = 2^(..) */
          ikk++; } } } /* kz,ky,kx */
    if (iQ!=Q->N || ikk!=Nkk) ERROR(("internal: Q->N=%d (%d expected)  Nkk=%d (%d expected)\n\
*** (may be caused by K rounding problems - try to change K a bit)",iQ,Q->N, ikk,Nkk))

    VV(oldL,=boxL) }

  /*** preparing tables of multiplying factors ex=exp[rx*..], ey, ez, etc.
       and the total dipole moment ***/

  if (Nq==0) return 0; /* nothing to do */

  qtabsize=Nq*sizeof(qtab_t);

  if (pass==1) {
    VO(Q->M,=0)
    size=Nq*sizeof(complex);
    /* ex,ey,ez will be freed in pass 2 */
    alloc(qtab,qtabsize);
    alloc(ex,size); alloc(ey,size); alloc(ez,size); }
  else {
    allocarrayzero(F,Nq); }

  /*** configuration -> Ewald interface ***/

  iq=0;
  if (pass==1)
    loop (n,0,el.No) {
//prt("%d %g\n",n,q[n]);
      if ((aux=q[n])) {
        VV(Q->M,+=aux*rp[n])
        qtab[iq].x.re=aux; qtab[iq].x.im=0;
        ex[iq].re=cos(aux=PIL[0]*rp[n][0]); ex[iq].im=sin(aux);
        ey[iq].re=cos(aux=PIL[1]*rp[n][1]); ey[iq].im=sin(aux);
        ez[iq].re=cos(aux=PIL[2]*rp[n][2]); ez[iq].im=sin(aux); iq++; }
    } /* n */
  else /* pass 2 */
    loop (n,0,el.No) {
      if ((aux=q[n])) {
        /* normal charge (Drude may follow) */
        qtab[iq].x.re=aux; qtab[iq].x.im=0; iq++; }
    }
    if (iq!=Nq) ERROR(("Ewald:iq=%d Nq=%d el.No=%d",iq,Nq,el.No))                 


/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %             main loop over k-vectors  (passes 1 and 2)              %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

  Qk=Q->Qk;
  iQ=ikk=0;
  kxmax=(int)KX;

  /* serial version */

  if (pass==1) {
#include "ewpass1.c"
  }
  if (pass==2) {
    if (el.measure) {
#if PRESSURETENSOR&PT_VIR
      double PVIR[PT_DIM];
      double ptf,E1;

      memset(PVIR,0,sizeof(PVIR));
#endif /*# PRESSURETENSOR&PT_VIR */
#include "ewpass2m.c"
      Energy/=4*PI;
#if PRESSURETENSOR&PT_VIR
      VO(el.Pvir,+=Energy)
      loop (kx,0,PT_DIM) el.Pvir[kx]-=2/(4*PI)*PVIR[kx];
#endif /*# PRESSURETENSOR&PT_VIR */
    } else {
#include "ewpass2.c"
    } }

  /* end of serial version */
  /* because the code below calculating the M^2-term constributions */
#define PVIR el.Pvir

  if (pass==2) {
    double Edip;

    VO(ek,=0) /* dipolar force */

    if (el.slab) {
      /* interface in slab pseudo-2D b.c. */
      /* see In-Chul Yeh and Max L. Berkowitz, JCP 111, 3155 (1999) */
      c = 4*PI/PROD(boxL);
      ek[2] = -c*Q->M[2];
      Edip = c/2*Sqr(Q->M[2]); }
    else {
      /* isotropic dielectric b.c. */
      c = 4*PI/(2*el.epsinf+1)/PROD(boxL);
      VV(ek,=-c*Q->M)
      Edip = c/2*SQR(Q->M); }

#if PRESSURETENSOR&PT_VIR
    VO(PVIR,+=Edip)
#endif /*# PRESSURETENSOR&PT_VIR */
    iq=0;

#if PRESSURETENSOR&PT_VIR
#  if PRESSURETENSOR&PT_OFF
#    define ADDFORCES(R) \
       VVV(ptr,+=F[iq],+aux*ek) \
       if (el.measure) { \
         VVV(PVIR,+=aux*ek,*R) \
         PVIR[3]+=aux*ek[1]*R[2]; \
         PVIR[4]+=aux*ek[2]*R[0]; \
         PVIR[5]+=aux*ek[0]*R[1]; }
#  else /*# PRESSURETENSOR&PT_OFF */
#    define ADDFORCES(R) \
       VVV(ptr,+=F[iq],+aux*ek) \
         if (el.measure) { \
           VVV(PVIR,+=aux*ek,*R) }
#  endif /*#!PRESSURETENSOR&PT_OFF */
#else /*# PRESSURETENSOR&PT_VIR */
#  define ADDFORCES(R) \
       VVV(ptr,+=F[iq],+aux*ek)
#endif /*#!PRESSURETENSOR&PT_VIR */

    if (!el.sf) {
      /* summing up forces and elst. potential, adding M-term */
      loop (n,0,el.No) {
#if PRESSURETENSOR&PT_VIR
#endif /*# PRESSURETENSOR&PT_VIR */
        if ((aux=q[n])) {
            ptr=frp[n];
            ADDFORCES(rp[n])
            iq++; }
      }
      if (iq!=Nq) ERROR(("Ewald internal iq=%d Nq=%d",iq,Nq)) }

    free(F); free(qtab);
    free(ez); free(ey); free(ex);

    if (el.kappa==0) return 0; /* patch !!! */

    if (el.measure) {
      double sumq=charges.sumq;

      return
        Energy
        - alsPI*sumq 
        + Edip; }
    else
      return 0;
  }

  return (double)Q->N/el.No;

} /* Ewald */

#undef CxE
#undef CxiE
#undef Cadd
#undef Csqr
#undef iqeloop
#undef iqfloop
#undef iqloop

/*
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                  Ewald tests and parameter setting                    %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

static double Ewaldparm(double cutoff,vector L) /***************** Ewaldparm */
{
  /* calculates alpha, kappa from cutoff, epsr, epsk */
  double x,kappa=el.kappa; /* ??? (for the case if K is single) */
  double V=PROD(L);
  int nit=10000;
/*.....  double minL=fmin(fmin(L[1],L[0]),L[2]);*/

  if (charges.sumq<=0) {
    /* just calculate charges */
    el.kappa=1; 
    Ewald_k_space(-1,NULL,L,NULL,NULL); 
    el.kappa=kappa; }

  if (charges.sumq>0) {
    /* alpha, random charges */
    el.alpha=log(2*charges.max*sqrt(charges.sumq/cutoff/V)/el.epsr);
    if (el.alpha<0) {
      ERROR(("Ewald auto set: cannot determine alpha (check epsr)\n\
*** (charges.max=%g charges.sumq=%g cutoff=%g V=%g el.epsr=%g)",charges.max,charges.sumq,cutoff,V,el.epsr));
      return -1; }
    el.alpha=sqrt(el.alpha)/cutoff;
    /* kappa, random charges */
    do {
      if (!nit--) {
        ERROR(("Ewald auto set: cannot determine kappa (check el.epsk)\n\
*** (charges.max=%g el.alpha=%g el.epsk=%g)",charges.max,el.alpha,el.epsk))
        return -1; }
      x=kappa;
      kappa=log(charges.max*el.alpha/(PI*el.epsk)*sqrt(8*charges.sumq/(kappa*V)));
      if (kappa<=0) ERROR(("Ewald auto set: cannot determine kappa (check epsk)"))
      kappa=el.alpha*sqrt(kappa)/PI;
    } while (fabs(x-kappa)>1e-9);
    el.kappa=kappa; }

  return 0;
}

                                                          /*** Ewald_r_space */
double Ewald_r_space(vector *frp,vector boxL,vector *rp,real *q)
/* r-space Ewald sum - simple nearest-image code
   - library erfc used
   - cutoff<=L/2 required
   - normalized configuration (in basic cell) required
   - optimized for cube and cutoff~L/2
   - always measure (as if el.measure=1)
*/
{
#define ER 1.1283791670955126 /* 2/\pi */
  int i,j,l;
  vector Lh;
  double Ereal=0,er;
  double cq=Sqr(el.cutoff);

  VV(Lh,=0.5*boxL)

  loop (l,0,3) if (el.cutoff>Lh[l]*1.00001) ERROR(("cutoff=%g>L%c/2=%g",el.cutoff,l+'x',Lh[l]))

  loop (i,0,el.No) loop (j,i+1,el.No) {
    double rr,r1;
    double f,er1;
    vector dr;

    dr[0]=rp[j][0]-rp[i][0]; if (dr[0]<-Lh[0]) dr[0]+=boxL[0]; else if (dr[0]>Lh[0]) dr[0]-=boxL[0];
    dr[1]=rp[j][1]-rp[i][1]; if (dr[1]<-Lh[1]) dr[1]+=boxL[1]; else if (dr[1]>Lh[1]) dr[1]-=boxL[1];
    dr[2]=rp[j][2]-rp[i][2]; if (dr[2]<-Lh[2]) dr[2]+=boxL[2]; else if (dr[2]>Lh[2]) dr[2]-=boxL[2];

    rr=SQR(dr);
    if (rr>=cq) continue;

    r1=sqrt(rr);
    er1=erfc(el.alpha*r1);
    Ereal+=q[i]*q[j]*er1/r1;
    f=q[i]*q[j]/(rr*r1)*(ER*(el.alpha*r1)*exp(-Sqr(el.alpha*r1))+er1);

    VV(frp[i],-=f*dr)
    VV(frp[j],+=f*dr)
#if PRESSURETENSOR&PT_VIR
    VVV(PVIR,+=f*dr,*dr)
#  if PRESSURETENSOR&PT_OFF
    PVIR[3]+=f*dr[1]*dr[2];
    PVIR[4]+=f*dr[2]*dr[0];
    PVIR[5]+=f*dr[0]*dr[1];
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_VIR */
  }
  return Ereal;
}

                                                   /*** Ewald_r_space_spline */
double Ewald_r_space_spline(vector *frp,vector boxL,vector *rp,real *q)
/* the same as above but using spline version (use initerfc first) */
{
#define ER 1.1283791670955126 /* 2/\pi */
  int i,j,l;
  vector Lh;
  double Ereal=0,er;
  ermacrodcl /* to use macros eru, erd, erud */
  double cq=Sqr(el.cutoff);

  VV(Lh,=0.5*boxL)

  loop (l,0,3) if (el.cutoff>Lh[l]*1.00001) ERROR(("cutoff=%g>L%c/2=%g",el.cutoff,l+'x',Lh[l]))

  loop (i,0,el.No) loop (j,i+1,el.No) {
    double rr,r1;
    double f,er1;
    vector dr;

    dr[0]=rp[j][0]-rp[i][0]; if (dr[0]<-Lh[0]) dr[0]+=boxL[0]; else if (dr[0]>Lh[0]) dr[0]-=boxL[0];
    dr[1]=rp[j][1]-rp[i][1]; if (dr[1]<-Lh[1]) dr[1]+=boxL[1]; else if (dr[1]>Lh[1]) dr[1]-=boxL[1];
    dr[2]=rp[j][2]-rp[i][2]; if (dr[2]<-Lh[2]) dr[2]+=boxL[2]; else if (dr[2]>Lh[2]) dr[2]-=boxL[2];

    rr=SQR(dr);
    if (rr>=cq) continue;

    er=erud(rr);
    Ereal+=q[i]*q[j]*er;
    f=q[i]*q[j]*byerd;

    VV(frp[i],-=f*dr)
    VV(frp[j],+=f*dr)
#if PRESSURETENSOR&PT_VIR
    VVV(PVIR,+=f*dr,*dr)
#  if PRESSURETENSOR&PT_OFF
    PVIR[3]+=f*dr[1]*dr[2];
    PVIR[4]+=f*dr[2]*dr[0];
    PVIR[5]+=f*dr[0]*dr[1];
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_VIR */
  }
  return Ereal;
}

/*****************************************************************************
  fast (but memory consuming) method to calculate the following functions:

  eru(x) = beta*alpha*e(sqrt(x)*alpha)
  erd(x) = -beta*alpha^3*e'(sqrt(x)*alpha)/(sqrt(x)*alpha)

where

  e(y) = erfc(y)/y, e'(y)=de(y)/dy,

                         infinity
  erfc(y) = 2/sqrt(PI) * integral exp(-t**2) dt,
                            y

and alpha and beta are scaling parameters

eru,erd are macros in this version: use erud and byerd for both

method:
  hyperbolic spline (functions A+B/(x+C), no jumps in the derivative)

initialization:
  initerfc(int Grid, double minr, double maxr, double cutoff, double alpha, double beta, int shift)
where:
  Grid is a number of grid points per unity; negative Grid = less verbose
  [minr,maxr) defines the range of r, the range of the argument x of eru, erd,
    and erud is then [minr^2,maxr^2)
  cutoff and shift: 
    shift&1: eru is shifted to avoid jump at r=cutoff
    shift&2: erd is shifted to avoid jump at r=cutoff
    shift&4: table of eru,erd is printed to ertest.dat
  cutoff<=rmax required, usually cutoff=rmax

note: r < minr leads to low-accuracy results
      r > maxr is out of table 

*****************************************************************************/

#define DIM 3

struct Erfc_s Erfc = { DIM,4 };

real byerd;

static void erfcs(long to,double *e,double *z) /********************** erfcs */
/*
  e=eru(h*to) and z=erd(h*to) are computed.
  Must be called with decreasing sequence of parameters to, starting
  with a sufficiently high value because the integral is for this to
  approximated by the asymptotic formula
*/
{
#define f(X) (sqrt(X)*exp(X))
#define SQRTPI 1.772453850905516
  static long t= -1;
  static double sum;
  double
    h=1e0/(Erfc.ndiv*Erfc.sgrid),
    w=h/2, ha=h*0.2113248654051871, hb=h-ha,
    x=to*h;

  if (x>706) ERROR(("erfcs: range too large - cannot calculate splines"))
  else if (x>67) {
    static int pass;
    if (pass==0) {
      WARNING(("erfcs: range unnecessarily large (er<1e-30)\n\
*** (more warnings suppresed)"))
      pass=1; } }

  if (to>t) {
    t=to; x=to*h;
    if (x<1) ERROR(("erfcs: range too short"))
    /* the correct asymptotic formula is
       exp(-x)/sqrt(x+1-0.75/x+2/x^2-7.75/x^3+-... */
    sum=exp(-x)/sqrt(x+1-0.7/x);
  }

  while (t>to) {
    x=h*t;
    sum += w/f(x-ha)+w/f(x-hb);
    t--; }


  x=t*h+1e-20;
  *e=sum/sqrt(x);
  *z=(*e+2/exp(x))/x/SQRTPI;
  *e /= SQRTPI;
#undef f
} /* erfcs */

static void eerfcs(double q, long ito) /***************************** eerfcs */
/* ito,ifrom in units of step */
{
  long x;
  double e,z,eh,zh,ee,zz,ed,zd,xh,xxh;
  double step=1e0/Erfc.sgrid, h=step/Erfc.ndiv;
  ertab_p p;

  x=ito*Erfc.ndiv; /* in units of h */
  erfcs(x+1,&eh,&zh);
  erfcs(x  ,&e, &z);
  erfcs(x-1,&ee,&zz);
  ed=(eh-ee)/(2*h)*q; zd=(zh-zz)/(2*h)*q;
  xh=x*h;

  while (ito) { ito--;
    xxh=xh; ee=e; zz=z;
    erfcs(x=ito*Erfc.ndiv,&e,&z);
    xh=x*h;

    p=Erfc.tab+ito;
    p->Cu=(xh*(ee-e)-xxh*ed*step) / (ed*step-(ee-e));
    p->Bu= -ed*(p->Cu+xxh)*(p->Cu+xxh);
    p->Au=ee+ed*(p->Cu+xxh);
    ed= -p->Bu/(p->Cu+xh)/(p->Cu+xh);

    p->Cd=(xh*(zz-z)-xxh*zd*step) / (zd*step-(zz-z));
    p->Bd= -zd*(p->Cd+xxh)*(p->Cd+xxh);
    p->Ad=zz+zd*(p->Cd+xxh);
    zd= -p->Bd/(p->Cd+xh)/(p->Cd+xh); }
} /* eerfcs */

void initerfc(int Grid, double minr, double maxr, double cutoff, double alpha, double beta, int shift)
/* Grid<0 means suppressed info */
{
  alpha=fabs(alpha); /* negative = verbose */
  {
    long ifrom,ito,i;
    int grid=abs(Grid),verbose=Grid>0;
    double x,e,ee,z,zz;
    double h=1e0/(Erfc.ndiv*grid);
    double mine=0,maxe=0,minz=0,maxz=0,minre=0,maxre=0,minrz=0,maxrz=0;
    double Ad=beta*Cub(alpha), AuBd=beta*alpha, Bu=beta/alpha, Cud=1/Sqr(alpha);
    ertab_p er_p;

    ifrom=(int)(Sqr(minr*alpha)*grid); ito=(int)(Sqr(maxr*alpha)*grid)+1;
    /* ifrom,ito are in units of 1/grid */

    if (Erfc.grid==grid && Erfc.shift==shift
        &&Erfc.minr==minr && Erfc.maxr>=maxr
        && Erfc.alpha==alpha && Erfc.beta==beta) {
      Erfc.to=Erfc.tab+ito;
      prt(":::::: erfc already initialized ::::::\n");
      return; }

    Erfc.shift=shift;
    Erfc.sgrid=(real)grid;
    Erfc.grid=grid; Erfc.minr=minr; Erfc.maxr=maxr;
    Erfc.alpha=alpha; Erfc.beta=beta;

    Erfc.A=beta*alpha*(-2/SQRTPI);
    Erfc.alphaq=Sqr(alpha);
    Erfc.B=Erfc.A*Erfc.alphaq*2;

    if (Erfc.tab!=NULL) free(Erfc.tab);
    allocarray(Erfc.tab,ito);

    prt("\n:::::: erfc :::::: tab 0..%ld..%ld = %ld B  grid=%d\n",
        ifrom,ito,ito*sizeof(ertab_t),grid);

    if (verbose) prt("subgrid=%d  alpha=%g  beta=%g  range=[%g,%g)  cutoff=%g\n",
                     Erfc.ndiv,alpha,beta,minr,maxr,cutoff);

    Erfc.from=Erfc.tab+ifrom; Erfc.to=Erfc.tab+ito;

    eerfcs(1.0,ito);

    /***** check accuracy *****/
    ifrom *= Erfc.ndiv; ito *= Erfc.ndiv;
    /* ifrom,ito are in units of h */
    for (i=ito-1; i>=ifrom; i--) {
      erfcs(i,&e,&z);
      x=i*h; ee=eru(x)-e; zz=erd(x)-z;
      if (i==ito-1)
        if (verbose) prt("eru(%f)=%e  erd(%f)=%e\n",x,e*AuBd,x,z*Ad);
      if (ee<mine) mine=ee; if (ee>maxe) maxe=ee;
      if (zz<minz) minz=zz; if (zz>maxz) maxz=zz;
      ee /= e; zz /= z;
      if (ee<minre) minre=ee; if (ee>maxre) maxre=ee;
      if (zz<minrz) minrz=zz; if (zz>maxrz) maxrz=zz; }

    if (verbose) {
      prt("abs err:  eru in (%.3e,%.3e)  erd in (%.3e,%.3e)\n",
          mine*AuBd,maxe*AuBd, minz*Ad,maxz*Ad);
      prt("rel err:  eru in (%.3e,%.3e)  erd in (%.3e,%.3e)\n",
          minre,maxre, minrz,maxrz); }

    /* rescaling */
    loop (er_p,Erfc.tab,Erfc.to) {
      er_p->Ad *= Ad;   er_p->Au *= AuBd;
      er_p->Bd *= AuBd; er_p->Bu *= Bu;
      er_p->Cd *= Cud;  er_p->Cu *= Cud; }

    Erfc.sgrid *= Sqr(alpha);

    z=Sqr(cutoff);

    if (shift&1)  {
      x=eru(z);
      loop (er_p,Erfc.tab,Erfc.to) er_p->Au -= x;
      prt("eru(r) shifted by %g to avoid jump in elst energy\n",x); }

    if (shift&2)  {
      x=erd(z);
      loop (er_p,Erfc.tab,Erfc.to) er_p->Ad -= x;
      prt("erd(r) shifted by %g to avoid jump in elst forces\n",x); }

    if (shift&4)
#  include "ertest.c"
  }
} /* initerfc */

/***
    exact but slow variant of erud w/o direct elst term:
    - value of eru minus direct elst term returned,
    - erd minus direct elst term passed,

    exacterud_sqrt(x) = beta*alpha*[e(sqrt(x)*alpha)-1/(sqrt(x)*alpha)]
    exacterud_sqrt_1(x) = beta*alpha*[e(sqrt(x)*alpha)-1/(sqrt(x)*alpha)-2/sqrt(pi)]

    erd = -beta*alpha^3*[e'(sqrt(x)*alpha)+1/(sqrt(x)*alpha)^2]/(sqrt(x)*alpha)

    see above for e, e'

    efficient for low xx only -- to be used as Ewald correction for
    intramolecular charges close together, e.g. if a dipole is
    approximated by two point charges very close together.  Then the
    number of terms in the series is low and exacterud works reasonably
    fast.
***/

double exacterud_sqrt(double xx,real *erd)
{
  int n=1;
  double su=0,sd=0,q,x;

  xx*=Erfc.alphaq;

  q=1;

  while (fabs(q)>erexact_eps) {
    /*.....  if (fabs(q)>1e5) { Error("exacterud"); }*/
    x=q/(2*n+1);
    sd+=x*n;
    su+=x;
    q*=-xx/++n; }

  *erd=Erfc.B*sd;

  return Erfc.A*(1-su*xx);
}

double exacterud_sqrt_1(double xx,real *erd)
{
  int n=1;
  double su=0,sd=0,q,x;

  xx*=Erfc.alphaq;

  q=1;

  while (fabs(q)>erexact_eps) {
    /*.....  if (fabs(q)>1e5) { Error("exacterud"); }*/
    x=q/(2*n+1);
    sd+=x*n;
    su+=x;
    q*=-xx/++n; }

  *erd=Erfc.B*sd;

  return Erfc.A*(-su*xx);
}
/* end of module */

// 2D arrays probably wont work see: http://www.astro.gla.ac.uk/~norman/star/sc13/sc13.htx/N-a2b5c4.html
//void ewald_(double xyz[],double fx[][],double fy[][],double fz[][],double q[],double *energy,double *Lx,double *Ly,double *Lz,double *cutin,double *alfain,double *kappain,double *epsinf,int *natom,int *iw)
void ewald_(double xyz[],double grad[],double q[],double *energy,double *Lx,double *Ly,double *Lz,double *cutin,double *alfain,double *kappain,double *epsinf,int *natom,int *iw)
{
  static vector boxL;
  static vector *r,*f;
  double Ek=0.0e0,Er=0.0e0;
  double t0,t1,t2;
  int i,ipom2;
  double  sq;
  double fac=3.1577322e5*1.8897; // autoKK*ang from modules.f90

  if (*iw==0) {
   el.verbose=1;
   el.measure=1;
   el.cutoff= *cutin;
   el.alpha= *alfain;
   el.kappa= *kappain;
   el.epsinf= *epsinf;
   el.No=*natom;
   el.grid=0;  //zatim bez splajnu
   boxL[0] = *Lx;
   boxL[1] = *Ly;
   boxL[2] = *Lz;
   sq=0;
 
// loop (i,0,el.No) prt("charge %i = %g \n",i,q[i]);
// loop (i,0,el.No) q[i]*=sqrt(167100.75);
   loop (i,0,el.No) sq+=q[i];
   printf("sum q=%g\n",sq);
   if (fabs(sq)>1e-9) Error("nonzero sum of charges");
   sq/=el.No;
// TODO: pokud se zmeni naboje tady,tak je musime taky updatovat v nabu
   loop (i,0,el.No) q[i]-=sq;

   Ewald_k_space(0,NULL,boxL,NULL,q);
   if (el.grid) initerfc(el.grid,0.5,el.cutoff,el.cutoff,el.alpha,1,3);
 
   allocarrayzero(f,el.No);
 
   allocarray(r,el.No);

 }
// -------------END OF INIT------------------------------------
 
//-------------------ACTUAL COMPUATATION---------------------------
  if(*iw>0) {
	  ipom2=0;
  for(i=0;i<el.No;i++) {
	  r[i][0]=xyz[ipom2];
	  r[i][1]=xyz[ipom2+1];
	  r[i][2]=xyz[ipom2+2];
	  f[i][0]=0.0e0;
	  f[i][1]=0.0e0;
	  f[i][2]=0.0e0;
	  ipom2=ipom2+3;
  }
  Ewald_k_space(1,f,boxL,r,q);
  Ek=Ewald_k_space(2,f,boxL,r,q);

 /* t1=(double)clock()/CLOCKS_PER_SEC;
  if (el.grid)
    Er=Ewald_r_space_spline(f,boxL,r,q);
  else
    Er=Ewald_r_space(f,boxL,r,q);
  t2=(double)clock()/CLOCKS_PER_SEC;
*/
  /*
  prt("\nc=%g α=%g κ=%g Ek=%g Er=%g E=%.12g\n",
      el.cutoff,el.alpha,el.kappa,Ek,Er,Ek+Er);
*/
  //  prt("k-space time = %g s  r-space time = %g s  sum = %g s\n",t1-t0,t2-t1,t2-t0);

  ipom2=0;
  for(i=0;i<el.No;i++) {
	 grad[ipom2]=f[i][0];  //TODO jednodussi bude predavat rovnou fx,fy,fz
	 grad[ipom2+1]=f[i][1];
	 grad[ipom2+2]=f[i][2];
	 ipom2+=3;
  }
 //fx[iw][i]+=f[i][0]/fac;  //TODO jednodussi bude predavat rovnou fx,fy,fz
 //fy[iw][i]+=f[i][1]/fac;
 //fz[iw][i]+=f[i][2]/fac;
/*
  prt("FORCES:\n");
  loop (i,0,el.No)
    prt("%17.10f %17.10f %17.10f\n",VARG(f[i]));
*/
  *energy=Ek+Er;
}

}
 /*
  if (narg<3) {
    fprintf(stderr,"Ewald summation - test of function call. Call by:\n\
  ewaldf cutoff alpha kappa epsinf [grid] < DATA\n\
PARMS:\n\
  cutoff: r-space cutoff\n\
  alpha:  separation parameter, in A^-1\n\
  kappa:  max k-vector divided by 2 PI (K/box), in A^-1\n\
  epsinf: dielectric constant of the surrounding continuum (dipoles only)\n\
  grid:   grid=0 (default): use library erfc\n\
          grid>0: use splines\n\
ENV:\n\
  UNITS: any value = cook internal units (charges in e, r in A, energy in K)\n\
         not defined = as is (CGS); pressure in Pa is meaningless\n\
DATA = file of:\n\
  NS\n\
  Lx Ly Lz\n\
  x y z q [repeated NS times]\n\
lengths are in A, charges in e\n\
Examples:\n\
  ewaldf 9 .7 1.2 1  < ewald.dat\n\
  ewaldf 10 .6 1.1 1  < ewald.dat\n\
");
    exit(0); }
*/
//  t0=(double)clock()/CLOCKS_PER_SEC;
   /*
   loop (int i,0,el.No) {
     if (!gets(line)) Error("no data");
     sscanf(line,"%lf%lf%lf%lf",r[i],r[i]+1,r[i]+2,&q[i]); }
 */
