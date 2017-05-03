//Interface to the ewald function from MACSIMUS package,
//released under GNU GPL licence and kindly modified by author Jiri Kolafa.
typedef double real;
#include "EWALD/include.h"
#include "EWALD/vector3d.h"
#include "EWALD/ewaldf.h"
static  vector boxL;
static  vector *r,*f;
static  real *q;
struct el_s el;

int ewald_init(double *Lx,double *Ly,double *Lz,double *cutin,double *alfain,double *kappain,int *natom3,int *ipom)
{
  int i;
  double sq;
  el.cutoff=*cutin;
  el.alpha=*alfain;
  el.kappa=*kappain;
  el.epsinf=3e33;
  el.No=*natom3;
  boxL[0] = *Lx;
  boxL[1] = *Ly;
  boxL[2] = *Lz;

//  if (narg>4) el.epsinf=atof(arg[4]);
//  if (narg>5) el.grid=atof(arg[5]);
  Ewald_k_space(0,NULL,boxL,NULL,q);
  if (el.grid) initerfc(el.grid,0.5,el.cutoff,el.cutoff,el.alpha,1,3);

  allocarrayzero(f,el.No);

  allocarray(r,el.No);
  allocarray(q,el.No);

  /*
  loop (int i,0,el.No) {
    if (!gets(line)) Error("no data");
    sscanf(line,"%lf%lf%lf%lf",r[i],r[i]+1,r[i]+2,&q[i]); }
*/
  sq=0;
  loop (i,0,el.No) sq+=q[i];
  printf("sum q=%g\n",sq);
  if (fabs(sq)>1e-9) Error("nonzero sum of charges");
  sq/=el.No;
  loop (i,0,el.No) q[i]-=sq;

  loop (i,0,el.No) q[i]*=sqrt(167100.75);
  el.verbose=1;
  el.measure=1;
}
int ewald(double xyz[],double grad[],double *energy,int *ipom)
{
  double Ek=0.0,Er=0.0;
  double t0,t1,t2;
  int i;
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

  /* double Ewald(int pass,vector boxL,real *q,vector *frp,vector *rp); */

//  t0=(double)clock()/CLOCKS_PER_SEC;
  Ewald_k_space(1,f,boxL,r,q);
  Ek=Ewald_k_space(2,f,boxL,r,q);
 /* t1=(double)clock()/CLOCKS_PER_SEC;
  if (el.grid)
    Er=Ewald_r_space_spline(f,boxL,r,q);
  else
    Er=Ewald_r_space(f,boxL,r,q);
  t2=(double)clock()/CLOCKS_PER_SEC;
*/
  prt("\nc=%g α=%g κ=%g Ek=%g Er=%g E=%.12g\n",
      el.cutoff,el.alpha,el.kappa,Ek,Er,Ek+Er);
//  prt("k-space time = %g s  r-space time = %g s  sum = %g s\n",t1-t0,t2-t1,t2-t0);

  prt("FORCES:\n");
  loop (i,0,el.No)
    prt("%17.10f %17.10f %17.10f\n",VARG(f[i]));

  *energy=Ek;

}
