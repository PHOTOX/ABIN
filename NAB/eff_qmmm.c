/***********************************************************************
                            MMEQMMM()
************************************************************************/

/*
 * Here is the mme function for 3D or 4D, depending upon the dim variable.
 *
 * Calling parameters are as follows:
 *
 * x - input: the atomic (x,y,z) coordinates
 * f - updated: the gradient vector
 * iter - the iteration counter, which if negative selects the following:
 *        -1 print detailed energy values
 *        -3 call nbond or egb to deallocate static arrays, then deallocate
 *           grad
 *        -4 force a nonbonded list update (option for xmin/lmod)
 *        -(any other value) normal execution, but no printout and no
 *                           nonbonded list update (-2 is used in xmin/lmod)
 */

static
REAL_T mme_qmmm(REAL_T * x, REAL_T * f, int *iter)
{
   REAL_T ebh, eba, eth, eta, eph, epa, enb, eel, enb14, eel14, ecn, edssp;
   REAL_T e_gb, esurf, evdwnp, frms, enbips, eelips;
   REAL_T e_pb, e_rism;
   REAL_T ene[20];
   REAL_T tmme1, t1, t2, treduce1;
   int i, j, k, goff, threadnum, numthreads, maxthreads;
   int iteration, mask, consumer, producer, numcopies;
   int dummy = 0;
   size_t n;

   static REAL_T *grad = NULL;
   static int *reqack = NULL;


   t1 = seconds();
   tmme1 = t1;
   threadnum = mytaskid;
   numthreads = numtasks;
   maxthreads = 1;
   n = (size_t) prmqm->Natom;

   maxthreads = 1;

   numcopies = 1;


   /*
    * Build the non-bonded pair list if it hasn't already been built;
    * rebuild it every nsnb iterations.  The non-bonded pair list
    * uses blocksize to group OpenMP thread to loop index, or MPI task
    * to loop index, mapping in the nblist and egb functions.  It is
    * global and fully populated for OpenMP, and local and partially
    * populated for MPI and SCALAPACK.
    *
    * If Generalized Born surface area calculations are selected,
    * build the non-polar pair list if it hasn't already been built;
    * rebuild it every nsnp iterations.  This pair list will be used for
    * non-polar first derivative calculations (for line minimization
    * within Newton-Raphson) and for non-polar egb interactions.  It
    * is global and fully populated for OpenMP, and local and fully
    * populated for MPI and SCALAPACK, which although wasteful of
    * memory, allows the full flexibility of loop index to MPI process
    * mapping that occurs with dynamic loop scheduling.
    */


   /* hcp does not use pairlists */
   if (!hcp)
   {

      if (nb_pairs_qm < 0 || (*iter > 0 && *iter % nsnb == 0)
                    ||  *iter == -4) {
         t2 = seconds();
         *tmmeOther += t2 - t1;
         t1 = t2;
         if( prmqm->IfBox ){
            nb_pairs_qm = nblist_box( x, lpairs_qm, pairlist_qm, cut );
         } else {
            nb_pairs_qm = nblist(lpairs_qm, upairs_qm, pairlist_qm, x, dummy, 1, cut,
                        prmqm->Natom, dim, frozen);
         }
         t2 = seconds();
         *tmmePair += t2 - t1;
         t1 = t2;
      }


   }


   if (grad == NULL) {
      grad = vector(0, numcopies * dim * n);
   }

   for (i = 0; i < dim * prmqm->Natom; i++) {
      grad[i] = 0.0;
   }

   t2 = seconds();
   *tmmeOther += t2 - t1;
   t1 = t2;

   ebh = ebond(prmqm->Nbonh, prmqm->BondHAt1, prmqm->BondHAt2,
               prmqm->BondHNum, prmqm->Rk, prmqm->Req, x, grad);
   eba = ebond(prmqm->Mbona, prmqm->BondAt1, prmqm->BondAt2,
               prmqm->BondNum, prmqm->Rk, prmqm->Req, x, grad);
   ene[3] = ebh + eba;
   t2 = seconds();
   *tmmeBond += t2 - t1;
   t1 = t2;

   eth = eangl(prmqm->Ntheth, prmqm->AngleHAt1, prmqm->AngleHAt2,
               prmqm->AngleHAt3, prmqm->AngleHNum, prmqm->Tk, prmqm->Teq, x, grad);
   eta =
       eangl(prmqm->Ntheta, prmqm->AngleAt1, prmqm->AngleAt2, prmqm->AngleAt3,
             prmqm->AngleNum, prmqm->Tk, prmqm->Teq, x, grad);
   ene[4] = eth + eta;
   t2 = seconds();
   *tmmeAngl += t2 - t1;
   t1 = t2;

   eph = ephi(prmqm->Nphih, prmqm->DihHAt1, prmqm->DihHAt2,
              prmqm->DihHAt3, prmqm->DihHAt4, prmqm->DihHNum,
              prmqm->Pk, prmqm->Pn, prmqm->Phase, x, grad);
   epa = ephi(prmqm->Mphia, prmqm->DihAt1, prmqm->DihAt2,
              prmqm->DihAt3, prmqm->DihAt4, prmqm->DihNum,
              prmqm->Pk, prmqm->Pn, prmqm->Phase, x, grad);
   ene[5] = eph + epa;
   ene[6] = 0.0;                /*  hbond term not in Amber-94 force field */
   t2 = seconds();
   *tmmePhi += t2 - t1;
   t1 = t2;

   /* In the following lpairs is a dummy argument that is not used. */

   nbond(lpairs_qm, prmqm->N14pairs, N14pearlist_qm, 1, x, grad, &enb14, &eel14,
         scnb, scee, 0);
   ene[7] = enb14;
   ene[8] = eel14;
   t2 = seconds();
   *tmmeNonb += t2 - t1;
   t1 = t2;

   if (e_debug) {
      EXPR("%9.3f", enb14);
      EXPR("%9.3f", eel14);
   }

   if (nconstrained) {
      ecn = econs(x, grad);
      t2 = seconds();
      *tmmeCons += t2 - t1;
      t1 = t2;
   } else {
      ecn = 0.0;
   }

   if (wdssp != 0.0) {
      edssp = fdssp(x, grad);
      t2 = seconds();
      *tmmeCons += t2 - t1;
      t1 = t2;
   } else {
      edssp = 0.0;
   }
   ene[9] = ecn + edssp;

   t2 = seconds();
   *tmmeOther += t2 - t1;
   t1 = t2;
   if (gb) {
      t1 = t2;
   } else if (pbsa && mytaskid==0) {
      t1 = t2;
   } else {
      enbips = eelips = 0.;
      nbond(lpairs_qm,upairs_qm, pairlist_qm, 0, x, grad, &enb, &eel, 1.0, 1.0, 0);
      ene[1] = enb + enbips;
      ene[2] = eel + eelips;
      ene[10] = 0.0;
      ene[11] = 0.0;
      ene[12] = 0.0;
      if (e_debug) {
         EXPR("%9.3f", enb);
         EXPR("%9.3f", eel);
      }

      t2 = seconds();
      *tmmeNonb += t2 - t1;
      t1 = t2;
   }


   /* Here is no reduction of the grad array.  Copy it to the f array. */

   for (i = 0; i < dim * prmqm->Natom; i++) {
     f[i] = grad[i];
   }


   /* Calculate the rms gradient. */

   frms = 0.0;
   for (i = 0; i < dim * prmqm->Natom; i++)
      frms += f[i] * f[i];
   frms = sqrt(frms / (dim * prmqm->Natom));

   /* Calculate the total energy. */

   ene[0] = 0.0;
   for (k = 1; k <= 12; k++) {
      ene[0] += ene[k];
   }

   /*
    * Print the energies and rms gradient but only for task zero,
    * and only for positive values of the iteration counter.
    */

   if (mytaskid == 0) {
      if (*iter > -1 && (*iter == 0 || *iter % ntpr == 0)) {
         fprintf(nabout,
                 "ff:%6d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2e\n",
                 *iter, ene[0], ene[3] + ene[4] + ene[5],
                 ene[1] + ene[7], ene[2] + ene[8],
                 ene[9] + ene[11] + ene[12], ene[10], frms);
         fflush(nabout);
      }

      /* A value of -1 for the iteration counter is reserved for printing.
       * Values of 2 and 3 for e_debug used for sander-like output printing 
       * for GB and PB, respectively: JMS 2/2011
       */

      if (*iter == -1) {
        fprintf(nabout, "     bond:  %15.9f\n", ene[3]);
        fprintf(nabout, "    angle:  %15.9f\n", ene[4]);
        fprintf(nabout, " dihedral:  %15.9f\n", ene[5]);
        fprintf(nabout, "    enb14:  %15.9f\n", ene[7]);
        fprintf(nabout, "    eel14:  %15.9f\n", ene[8]);
        fprintf(nabout, "      enb:  %15.9f\n", ene[1]);
        fprintf(nabout, "      eel:  %15.9f\n", ene[2]);
        fprintf(nabout, "      egb:  %15.9f\n", ene[10]);
        fprintf(nabout, "    econs:  %15.9f\n", ene[9]);
        fprintf(nabout, "    esurf:  %15.9f\n", ene[11]);
        fprintf(nabout, "    Total:  %15.9f\n", ene[0]);
        fprintf(nabout, "     frms:  %15.9f\n", frms);
      } 
      if (e_debug == 2) {
        fprintf(nabout, " BOND    = %13.4lf  ANGLE   = %13.4lf  DIHED      = %13.4lf\n", ene[3], ene[4], ene[5]);
        fprintf(nabout, " VDWAALS = %13.4lf  EEL     = %13.4lf  EGB        = %13.4lf\n", ene[1], ene[2], ene[10]);
        fprintf(nabout, " 1-4 VDW = %13.4lf  1-4 EEL = %13.4lf  RESTRAINT  = %13.4lf\n", ene[7], ene[8], ene[9]);
        fprintf(nabout, " ESURF   = %13.4lf\n", ene[11]);
    } else if (e_debug == 3) {
        fprintf(nabout, " BOND    = %13.4lf  ANGLE   = %13.4lf  DIHED      = %13.4lf\n", ene[3], ene[4], ene[5]);
        fprintf(nabout, " VDWAALS = %13.4lf  EEL     = %13.4lf  EPB        = %13.4lf\n", ene[1], ene[2], ene[10]);
        fprintf(nabout, " 1-4 VDW = %13.4lf  1-4 EEL = %13.4lf  RESTRAINT  = %13.4lf\n", ene[7], ene[8], ene[9]);
        fprintf(nabout, " ECAVITY = %13.4lf  EDISPER = %13.4lf", ene[11], 0.0);
      }
   }

   /* If static_arrays is 0, deallocate the gradient and reqack arrays. */

   if (!static_arrays) {
     if (grad != NULL) {
       free_vector(grad, 0, numcopies * dim * n);
       grad = NULL;
     }
     if (reqack != NULL) {
       free_ivector(reqack, 0, numcopies);
       reqack = NULL;
     }
   }

   /* Update the mme time, which includes the reduction time. */

   t2 = seconds();
   *tmmeOther += t2 - t1;
   *tmme += t2 - tmme1;
   return (ene[0]);
}
