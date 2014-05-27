#  define PT_VIR 1 /* virial (=configurational + constraint) part */
#  define PT_OFF 4 /* also off-diagonal terms (with PT_VIR and PT_KIN) */

#    if PRESSURETENSOR&PT_OFF
#      define PT_DIM (2*DIM)
#    else /*# PRESSURETENSOR&PT_OFF */
#      define PT_DIM (DIM)
#    endif /*#!PRESSURETENSOR&PT_OFF */

/*** sums of charges etc. ***/
extern struct charges_s {
  real sum,max,sumq;
} charges;

typedef struct { real re,im; } complex;

typedef struct /*SDS*/ {
  int size,N;
  vector M;         /* dipole moment */
  complex Qk[1];    /* array[Q_t.N] : sums Q */
} Q_t;

extern Q_t *Q;      /* Ewald pass 1 -> pass 2 interface variable */

/*** moved here in this interface ***/
extern struct el_s {
#if PRESSURETENSOR
  real Pvir[PT_DIM]; /* pressure tensor:
                        Pvir[0]=Pxx, Pvir[1]=Pyy, Pvir[2]=Pzz,
                        Pvir[3]=Pyz, Pvir[4]=Pzx, Pvir[5]=Pzy */
#endif

  real alpha;      /* Ewald separation parameter in units of 1/A;
                      also cutoff electrostatic parameter */
  real kappa;      /* k-space cutoff in form kappa=K[i]/L[i] */
  real cutoff;     /* r-space cutoff */
  real epsr,epsk;  /* expected r- and k-space errors in forces, in K/A */
                   /* not active here */
  real minqq;      /* minimum expected charge-charge separation
                      (to test Drude and to report error of the erfc approximation) */
  real epsinf;     /* dielectric constant at infinity
                      (must be very big if free charges are present */
  real diff;       /* max rel.change of L before a warning about "unchanged K" */
  real Perr;       /* pressure tensor threshold (basic; WARNING for 10*Perr) */
  real rplus;      /* additional range to erfc splines (former ERFCPLUS)
                      needed for optimized water models, 
                      does not hurt if unnecessarily large */
  //  real sat;        /* target saturation (see tau.sat) */
  int measure;     /* measurements on / off (a bit faster) */
  int No;          /* number of atoms */
  int verbose;     /* more prints */
  int rshift;      /* &1: eru(r)=erfc(r)/r spline shifted to avoid jump at cutoff
                      &2: as above for the derivative, erd(r) 
                      &4: prints eru(r), erd(r), d eru(r)/dr verbosely on init */
  int grid;        /* # of grid points per unity of squared distance */
  //  int test;        /* test module switch (former etest) */
  //  int diag;        /* whether to include diagonal correction */
  //  int centroid;    /* 1=dipole moments (of ions) are w.r.t charge centroid */
  int sf;          /* structure factor switch:
                      0 = standard k-space Ewald sum
                      1 = sphericalized (cube only)
                      3 = 3D */
  int slab;        /* slab correction, see
                      In-Chul Yeh and Max L. Berkowitz, JCP 111, 3155 (1999) */

  /* private: */
  int sfsize;      /* to check overflow */
  struct sfr_s {   /* for sphericalized (radial) structure factor */
    int nk;           /* # of independent vectors */
    real q;           /* squared structure factor */
  } *sfr;
  struct sf3d_s {  /* for 3D structure factor */
    int k[DIM];       /* k-vector */
    real q;           /* squared structure factor */
  } *sf3d;
  double sumM;     /* sum of molecule dipole moments = max.cell moment
                      for saturation autoset, one cfg value (cf. Eext.sumM)
                      POLAR: includes induced moments */
  //  double xinf;     /* 1/(el.epsinf*2+1) */
} el;


/*** erfc splines, macro version ***/

typedef struct { real Au,Bu,Cu,Ad,Bd,Cd; } ertab_t,*ertab_p;

#define erexact_eps 1e-10

extern struct Erfc_s {
  int key;          /* only 3 in this version */
  int ndiv;         /* to def integration step; must be at least 4 */
  int grid;         /* grid points/1 */
  int shift;        /* shiff to be continuous (no even tiny error) */
  ertab_p tab;      /* table offset */
  ertab_p from,to;  /* valid range of the table */
  /* private: */
  real sgrid;
  real alphaq,A,B;
  real minr,maxr,alpha,beta;
} Erfc;

extern real byerd; /* by-product result of macro erud
                      nb: local byerd may be faster, cf. ermacrodcl */

#define eru(x) (er_p=Erfc.tab+(int)((x)*Erfc.sgrid), \
  er_p->Bu/(er_p->Cu+(x))+er_p->Au)
#define erd(x) (er_p=Erfc.tab+(int)((x)*Erfc.sgrid), \
  er_p->Bd/(er_p->Cd+(x))+er_p->Ad)
#define erud(x) (er_p=Erfc.tab+(int)((x)*Erfc.sgrid), \
  byerd = er_p->Bd/(er_p->Cd+(x)) + er_p->Ad,  \
  er_p->Bu/(er_p->Cu+(x)) + er_p->Au)

#define ermacrodcl ertab_p er_p; real byerd;

void initerfc(int Grid, double minr, double maxr, double cutoff, double alpha, double beta, int shift);

double exacterud_sqrt(double xx,real *erd);
