/*
   Ewald pass 2 (calculate forces, energy, pressure tensor)
   #included from ewald.c

   Overview of arrays:
   ^^^^^^^^^^^^^^^^^^^
   ekk[ikk] = 4 exp[-pi^2 alpha^2 (k/L)^2] / [V (k/L)^2]
     ikk runs over k-vectors k=(kx,ky,kz), in the order defined by the 
     kx,ky,kz loops with all +-kx,+-ky,+-kz stored in the same index
     NB: factor 4 above comes from: 2 from differentiation of Q^2
                                    2 from counting k and -k at once 
   Qk[iQ] = sum_j q_j exp(2 pi i (k/L) . r_j]
     iQ runs over k-vectors k=(kx,ky,kz), in the order defined by the 
     kx,ky,kz loops; iQ runs faster than ikk
   qtab[iqt] depends on POLAR&32:
     iqt runs over [a subset of if PARALLEL] charges
     !(POLAR&32): qtab[iqt].x = q[iqt]*exp[2*PI*i*x*kx/Lx], etc;
     POLAR&32:    qtab[iqt].q = q[iqt]
                  qtab[iqt].x = exp[2*PI*i*x*kx/Lx], etc. (see ewald.c)

   NB: Energy, PVIR, Qk are in the PARALLEL=1|3 version macros
   referring to thread-based instances
*/

  jkk=0; /* for 3D structure factor only */

  loopto (kx,0,kxmax) {
    /* NOTE:
       calculations of qtab are unnecessarily repeated in the 2nd pass,
       however, storing all the factors would require a lot of memory,
       namely [ikk] indexed array of qtab's */
    if (kx) iqeloop(ex) CxE(qtab[iqt].x)
    iqloop qtab[iqt].xy = qtab[iqt].xY = qtab[iqt].x;

    /* kymax=(int)(L[1]*sqrt(Kk-Sqr(kx/L[0]))); */
    kymax=ktab[kx].kymax;

    loopto (ky,0,kymax) {
      if (ky) iqeloop(ey) { CxE(qtab[iqt].xy) CxiE(qtab[iqt].xY) }

      iqloop {
        qtab[iqt].xyz = qtab[iqt].xyZ = qtab[iqt].xy;
        qtab[iqt].xYz = qtab[iqt].xYZ = qtab[iqt].xY; }

      /* kzmax=(int)(L[2]*sqrt(Kk-Sqr(ky/L[1])-Sqr(kx/L[0]))); */
      kzmax=ktab[kx].kzmax[ky];
      loopto (kz,!(kx+ky),kzmax) {
        if (kz) iqeloop(ez) {
          CxE(qtab[iqt].xyz) CxiE(qtab[iqt].xyZ)
          CxE(qtab[iqt].xYz) CxiE(qtab[iqt].xYZ) }

        if (el.sf==0) {
          /********************************/
          /**    pass 2 - normal Ewald   **/
          /********************************/
          e=ekk[ikk];

#if PRESSURETENSOR&PT_VIR
/* to make more tables?...*/
	  ek[0]=e*(kx/boxL[0]); ek[1]=e*(ky/boxL[1]); ek[2]=e*(kz/boxL[2]);
	  kk=Sqr(kx/boxL[0])+Sqr(ky/boxL[1])+Sqr(kz/boxL[2]);
	  ptf=1/kk-PIaq;
#define CALCULATEPRESSURETENSOR \
          PVIR[0]+=E1*ptf*Sqr(kx/boxL[0]); \
          PVIR[1]+=E1*ptf*Sqr(ky/boxL[1]); \
          PVIR[2]+=E1*ptf*Sqr(kz/boxL[2]);
#define ASSIGN_E1 E1=
#else
#define CALCULATEPRESSURETENSOR /* empty */
#define ASSIGN_E1 /* empty */
	  ek[0]=e*kx/boxL[0]; ek[1]=e*ky/boxL[1]; ek[2]=e*kz/boxL[2];
#endif

          switch ((kx>0) + (ky>0) + (kz>0)) {

/* +++ */   case 3:
              A=Qk[iQ++]; B=Qk[iQ++]; C=Qk[iQ++]; D=Qk[iQ++];
              Energy += ASSIGN_E1 (Csqr(A)+Csqr(B)+Csqr(C)+Csqr(D))*e;
              CALCULATEPRESSURETENSOR
#if PRESSURETENSOR&PT_OFF
              PVIR[3]+=e*ptf*(Csqr(A)-Csqr(B)-Csqr(C)+Csqr(D))*(ky/boxL[1])*(kz/boxL[2]);
              PVIR[4]+=e*ptf*(Csqr(A)-Csqr(B)+Csqr(C)-Csqr(D))*(kx/boxL[0])*(kz/boxL[2]);
              PVIR[5]+=e*ptf*(Csqr(A)+Csqr(B)-Csqr(C)-Csqr(D))*(kx/boxL[0])*(ky/boxL[1]);
#endif
              iqfloop {
                FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re; /* +++ */
                FB = B.re*qtab[iqt].xyZ.im - B.im*qtab[iqt].xyZ.re; /* ++- */
                FC = C.re*qtab[iqt].xYz.im - C.im*qtab[iqt].xYz.re; /* +-+ */
                FD = D.re*qtab[iqt].xYZ.im - D.im*qtab[iqt].xYZ.re; /* +-- */
#if defined(POLAR) && POLAR&32
                Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                             + B.re*qtab[iqt].xyZ.re + B.im*qtab[iqt].xyZ.im
                             + C.re*qtab[iqt].xYz.re + C.im*qtab[iqt].xYz.im
                             + D.re*qtab[iqt].xYZ.re + D.im*qtab[iqt].xYZ.im);
#endif
                F[iqt][0] += QTIMES ek[0]*((FA+FB)+(FC+FD));
                F[iqt][1] += QTIMES ek[1]*((FA+FB)-(FC+FD));
                F[iqt][2] += QTIMES ek[2]*(FA-FB+FC-FD); }
              break;

/* 0++ */
/* +0+ */
/* ++0 */   case 2:
              A=Qk[iQ++]; B=Qk[iQ++];
              Energy += ASSIGN_E1 (Csqr(A)+Csqr(B))*e;
              CALCULATEPRESSURETENSOR
/* ..+ */     if (kz) {
#if PRESSURETENSOR&PT_OFF
                PVIR[3]+=e*ptf*(Csqr(A)-Csqr(B))*(ky/boxL[1])*(kz/boxL[2]);
                PVIR[4]+=e*ptf*(Csqr(A)-Csqr(B))*(kx/boxL[0])*(kz/boxL[2]);
                PVIR[5]+=e*ptf*(Csqr(A)+Csqr(B))*(kx/boxL[0])*(ky/boxL[1]);
#endif
                iqfloop {
                  FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re; /* ..+ */
                  FB = B.re*qtab[iqt].xyZ.im - B.im*qtab[iqt].xyZ.re; /* ..- */
#if defined(POLAR) && POLAR&32
                  Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                               + B.re*qtab[iqt].xyZ.re + B.im*qtab[iqt].xyZ.im);
#endif
                  if (kx) F[iqt][0] += QTIMES ek[0]*(FA+FB);
                  else    F[iqt][1] += QTIMES ek[1]*(FA+FB);
                          F[iqt][2] += QTIMES ek[2]*(FA-FB); } }
/* ++0 */     else {
#if PRESSURETENSOR&PT_OFF
                PVIR[3]+=e*ptf*(Csqr(A)-Csqr(B))*(ky/boxL[1])*(kz/boxL[2]);
                PVIR[4]+=e*ptf*(Csqr(A)+Csqr(B))*(kx/boxL[0])*(kz/boxL[2]);
                PVIR[5]+=e*ptf*(Csqr(A)-Csqr(B))*(kx/boxL[0])*(ky/boxL[1]);
#endif
                iqfloop {
                  FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re; /* .+0 */
                  FB = B.re*qtab[iqt].xYz.im - B.im*qtab[iqt].xYz.re; /* .-0 */
#if defined(POLAR) && POLAR&32
                  Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                               + B.re*qtab[iqt].xYz.re + B.im*qtab[iqt].xYz.im);
#endif
                  F[iqt][0] += QTIMES ek[0]*(FA+FB);
                  F[iqt][1] += QTIMES ek[1]*(FA-FB); } }
              break;

/* +00 */
/* 0+0 */
/* 00+ */   case 1:
              A=Qk[iQ++];
              Energy += ASSIGN_E1 Csqr(A)*e;
              CALCULATEPRESSURETENSOR
#if PRESSURETENSOR&PT_OFF
              PVIR[3]+=e*ptf*Csqr(A)*(ky/boxL[1])*(kz/boxL[2]);
              PVIR[4]+=e*ptf*Csqr(A)*(kx/boxL[0])*(kz/boxL[2]);
              PVIR[5]+=e*ptf*Csqr(A)*(kx/boxL[0])*(ky/boxL[1]);
#endif
              iqfloop {
                FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re;
#if defined(POLAR) && POLAR&32
                Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im);
#endif
                if (kx)      F[iqt][0] += QTIMES ek[0]*FA;
                else if (ky) F[iqt][1] += QTIMES ek[1]*FA;
                else         F[iqt][2] += QTIMES ek[2]*FA; }
              break;

/* ??? */   default: ERROR(("Ewald switch (k=%d %d %d)",kx,ky,kz)) } }

        else if (el.sf==1) {
          /***********************************************/
          /**  pass 2 - sphericalized structure factor  **/
          /***********************************************/

          if (!iscube(boxL))
	    ERROR(("sphericalized version of the structure factor is implemented for a cube only"))

          jkk=kx*kx+ky*ky+kz*kz;
          if (jkk>el.sfsize) ERROR(("el.sfr[%d] out of range [%d]",jkk,el.sfsize))

          switch ((kx>0) + (ky>0) + (kz>0)) {

/* +++ */   case 3:
              A=Qk[iQ++]; B=Qk[iQ++]; C=Qk[iQ++]; D=Qk[iQ++];
              el.sfr[jkk].q += Csqr(A)+Csqr(B)+Csqr(C)+Csqr(D);
              el.sfr[jkk].nk += 4;
              break;

/* 0++ */
/* +0+ */
/* ++0 */   case 2:
              A=Qk[iQ++]; B=Qk[iQ++];
              el.sfr[jkk].q += Csqr(A)+Csqr(B);
              el.sfr[jkk].nk += 2;
              break;

/* +00 */
/* 0+0 */
/* 00+ */   case 1:
              A=Qk[iQ++];
              el.sfr[jkk].q += Csqr(A);
              el.sfr[jkk].nk++;
              break;

/* ??? */   default: ERROR(("Ewald switch (k=%d %d %d)",kx,ky,kz)) } }

        else if (el.sf==3) {
          /************************************/
          /**  pass 2 - 3D structure factor  **/
          /************************************/

          if (jkk>el.sfsize) ERROR(("el.sf3d[%d] out of range [%d]",jkk,el.sfsize))

          switch ((kx>0) + (ky>0) + (kz>0)) {

/* +++ */   case 3:
              el.sf3d[jkk].k[0]=kx;
              el.sf3d[jkk].k[1]=ky;
              el.sf3d[jkk].k[2]=kz;
              A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A);

              el.sf3d[jkk].k[0]=kx;
              el.sf3d[jkk].k[1]=ky;
              el.sf3d[jkk].k[2]=-kz;
              A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A);

              el.sf3d[jkk].k[0]=kx;
              el.sf3d[jkk].k[1]=-ky;
              el.sf3d[jkk].k[2]=kz;
              A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A);

              el.sf3d[jkk].k[0]=kx;
              el.sf3d[jkk].k[1]=-ky;
              el.sf3d[jkk].k[2]=-kz;
              A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A);
              break;

/* 0++ */
/* +0+ */
/* ++0 */   case 2:
              el.sf3d[jkk].k[0]=kx;
              el.sf3d[jkk].k[1]=ky;
              el.sf3d[jkk].k[2]=kz;
              A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A);

              if (kz) {
                el.sf3d[jkk].k[0]=kx;
                el.sf3d[jkk].k[1]=ky;
                el.sf3d[jkk].k[2]=-kz;
                A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A); }
              else {
                el.sf3d[jkk].k[0]=kx;
                el.sf3d[jkk].k[1]=-ky;
                el.sf3d[jkk].k[2]=kz;
                A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A); }
              break;

/* +00 */
/* 0+0 */
/* 00+ */   case 1:
              el.sf3d[jkk].k[0]=kx;
              el.sf3d[jkk].k[1]=ky;
              el.sf3d[jkk].k[2]=kz;
              A=Qk[iQ++]; el.sf3d[jkk++].q += Csqr(A);

              break;

/* ??? */   default: ERROR(("Ewald switch (k=%d %d %d)",kx,ky,kz)) } 
        }
        ikk++; } } } /* kz,ky,kx */

  if (iQ!=Q->N) ERROR(("Ewald: iQ=%d != Q->N=%d; L=(%g %g %g)",
                               iQ,      Q->N,    VARG(boxL)))
