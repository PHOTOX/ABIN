#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;

static MOLECULE_T *my;

static REAL_T x[2];

int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __it0001__;
static STRING_T *__st0001__ = NULL;
static STRING_T *__st0002__ = NULL;
static STRING_T *__st0003__ = NULL;
my = getpdb_prm( STEMP( __st0001__, "input.pdb" ), STEMP( __st0002__, "leapre.ff99" ), STEMP( __st0003__, "" ), ITEMP( __it0001__, 0 ) );
my = getpdb( "input.pdb", NULL );

readparm( my, "input.top" );

mme_init( my, NULL, "::Z", x, NULL );



	exit( 0 );
}
