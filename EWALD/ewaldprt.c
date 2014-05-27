/* ewaldprt.c is include file for ewald.c */

if (charges.sumq!=0 && el.kappa>0) {
  double V=PROD(boxL);
  double acq=Sqr(el.alpha*el.cutoff);
  double dEr=charges.sumq*sqrt(el.cutoff/2/V)/acq/exp(acq);
  double dfr=2*sqrt(charges.sumq/el.cutoff/V)/exp(acq);
  double eaL=exp(-Sqr(PI/el.alpha*el.kappa));
  double dEk=charges.sumq*el.alpha/(PI*PI*el.kappa*sqrt(V*el.kappa))*eaL;
  double dfk=sqrt(8*charges.sumq/el.kappa/V)*el.alpha/PI*eaL;
  double rat=1.05;

// TO BE RECONSIDERED
#define electron 408.779582 /* charge of electron in program units */
  if (!el.sf && el.alpha>=0 && Sqr(charges.sum)/charges.sumq>2e-8)
    ERROR(("system is charged: SUM q=%g e, sum q^2=%g e",
           charges.sum/electron,charges.sumq/Sqr(electron)))

  if (el.epsinf<0) WARNING(("\
The dielectric constant of the surrounding continuum is negative:\n\
***   el.epsinf=%g",el.epsinf))

  if (fabs(el.epsinf+0.5)<1e-5)
    ERROR(("el.epsinf=%g out of range (too large dipolar response)",el.epsinf))

  prt("\n++++++ Ewald initialized ++++++  alpha = %.6f/A   kappa = %.5f/A\n",
      el.alpha,el.kappa);

  if (el.verbose && el.alpha>0) { /* <<<<<<<<<<<<<<<<<<<<<<<<<<<< */

    prt(
        "Q->N = #{k,|k|<=K,k=-k} = %i   "
        "Nkk = #{k,|k|<=K,kx>=0,ky>=0,kz>=0} = %i\n",
        iQ,ikk);
    prt("# of q = %d  SUM(q^2) = %g K*A = %g e^2\n\
SUM(q) = %g\n",
        Nq,charges.sumq,charges.sumq/Sqr(electron), charges.sum);

    /*** error estimates ***/

    /* standard deviations of expected errors caused by cutoffs */
    prt("\n\
The following estimates are based on the assumption of random distribution of\n\
charges. They work quite well for ionic systems. For dipolar systems without\n\
free charges (e.g., water), the k-space estimates are too pessimistic and\n\
actual errors are roughly 10x smaller; thus, epsk=10*epsr is recommended.\n");
    /* diagonal (always positive) terms omitted by |k|<K cutoff */
    prt("* expected non-diagonal cutoff error of total energy:\n\
  r: %.2e K        k: %.2e K\n", dEr,dEk);
    prt("* expected cutoff error of force on maximum charge %g = %g e:\n\
  r: %.2e K/A      k: %.2e K/A\n",
        charges.max, charges.max/electron,
        charges.max*dfr, charges.max*dfk);

    // DHmod: currently el.epsr and el.epsk =0, which results in floating point exception
    //rat=charges.max*dfr/el.epsr;
    if (rat>1.03)
      prt("WARNING: estimated r-space error is %.3fx larger than limit epsr=%g K/A\n",rat,el.epsr);
  //  rat=charges.max*dfk/el.epsk;
    if (rat>1.03)
      prt("WARNING: estimated k-space error is %.3fx larger than limit epsk=%g K/A\n",rat,el.epsk);

    /* square roots of sums of squared errors over all charges */
    e=sqrt(charges.sumq);
    prt("* expected summarized standard cutoff error of forces:\n\
  r: %.2e K/A      k: %.2e K/A\n", estimr=e*dfr,estimk=e*dfk); // Bug fixed by TT

    prt("NB: summarized standard errors are defined as sqrt[SUM_i(f_i-f_iref)^2]\n\
    - there is no division by N because of problematic definition in mixtures\n");

  } /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

  loop (i,0,3)
    if (el.cutoff>boxL[i]*0.50000001)
      prt("WARNING cutoff/boxL[%d] = %.7f > 1/2\n",i,el.cutoff/boxL[i]); }
