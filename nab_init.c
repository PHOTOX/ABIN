//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
#include "nabtypes.h"

static int mytaskid, numtasks;
static MOLECULE_T *m,*qm;
static STRING_T *in1 = "input.pdb";

extern REAL_T alpha_pme,cut_pme;
extern INT_T ipbc,nsnb_pme;
extern INT_T ips_pme;

void nab_init_(REAL_T *alphain,REAL_T *cutoffin,INT_T *nsnbin,INT_T *pbcin,INT_T *ipsin)
{
   nabout = stdout; /*default*/
   mytaskid = 0; numtasks = 1;
   alpha_pme = *alphain;
   cut_pme = *cutoffin;
   nsnb_pme = *nsnbin;
   ipbc = *pbcin;
   ips_pme = *ipsin;

   static REAL_T dummy[2];
   m = getpdb("input.pdb", NULL );

   readparm(m, "input.top");

   mme_init(m, NULL, "::Z", dummy, NULL);
}
