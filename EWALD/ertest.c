/* DEBUG: table of eru, erd, and num.deriv of eru */
{
  double x,e,u1,u2,d=1e-4;
  FILE *f=fopen("ertest.dat","wt");

  fprintf(f,"#  grid=%d  alpha=%g  range=%g  cutoff=%g\n",
          Erfc.grid,Erfc.alpha,maxr,cutoff);
  fprintf(f,"#  r       eru         erd   deriv(eru)\n");

  for (x=1; x<maxr; x+=0.001) {
    u1=eru(Sqr(x+d));
    u2=eru(Sqr(x-d));
    e=erud(x*x);
    fprintf(f,"%.3f %.15g %.15g %.15g\n",x,e,byerd,(u2-u1)/(2*d*x) ); }

  fclose(f);
}
