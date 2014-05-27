// DHmod: added alpha_pme, to be used in nbond and nbond2 routines
// should be somehow initialized by mme_init or perhaps directly in nabinit.c

//DHmod
REAL_T alpha_pme=1.0;

// variables mainly for nblist to know if we habe PBC
// becouse nblist does not have direct access to topology
// ipbc=1 - PERIODIC BOUNDARY CONDITIONS
INT_T ipbc;
// BOX sizes
REAL_T Lx_box,Ly_box,Lz_box;
REAL_T L2x,L2y,L2z;

// nonbonded cutoff, probably stupid to put it here but..
// reason: so that we can set it in nab_init, similarly to alpha_pme
// it is later copied to cut in sff_mme_init routine
REAL_T cut_pme;
INT_T nsnb_pme;
// ISOTROPIC PERIODIC SUM 1 - both LJ and Coulomb, 2 - Coulomb, 3 - LJ
INT_T ips_pme;
