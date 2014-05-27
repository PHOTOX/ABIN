#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
//#include "pme.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;
static MOLECULE_T *m,*qm;
static STRING_T *pdbin = "input.pdb";
static STRING_T *in1 = "input.pdb";
static STRING_T *topin = "input.top";
static STRING_T *leapin = "leaprc.ff99SB";
static STRING_T *dummyin = "";

extern REAL_T alpha_pme,cut_pme;
extern INT_T ipbc,nsnb_pme;
extern REAL_T Lx_box,Ly_box,Lz_box;
extern INT_T ips_pme;

void nab_init_(REAL_T *alphain,REAL_T *cutoffin,INT_T *nsnbin,INT_T *pbcin,INT_T *ipsin, INT_T *iqmmm)
{

	nabout = stdout; /*default*/
	mytaskid=0; numtasks=1;
	alpha_pme = *alphain;
	cut_pme = *cutoffin;
	nsnb_pme = *nsnbin;
	ipbc = *pbcin;
	ips_pme = *ipsin;

static REAL_T dummy[2];
static INT_T __it0001__;
static STRING_T *__st0001__ = NULL;
static STRING_T *__st0002__ = NULL;
static STRING_T *__st0003__ = NULL;
//my = getpdb_prm( STEMP( __st0001__, "input.pdb" ), STEMP( __st0002__, "leaprc.ff99" ), STEMP( __st0003__, "" ), ITEMP( __it0001__, 0 ) );
m = getpdb( "input.pdb", NULL );

readparm( m, "input.top" );

mme_init( m, NULL, "::Z", dummy, NULL );

if(*iqmmm == 1) {
 	qm = getpdb( "input_qm.pdb", NULL );
	readparm( qm, "input_qm.top" );
	mme_init( qm, NULL, "::Z", dummy, NULL );
}

//debugDEBUG
/*FILE_T *file1;
file1 = fopen( "grad237.dat", "w" );
fclose(file1);*/
}

void nab_getbox_(REAL_T *x,REAL_T *y,REAL_T *z) {
	*x = Lx_box;
	*y = Ly_box;
        *z = Lz_box;
}	

