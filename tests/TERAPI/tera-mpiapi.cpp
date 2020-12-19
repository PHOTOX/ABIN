/* int.cpp */
/***********************************************************************/
/*                                                                     */
/*  This file is part of the TeraChem software suite                   */
/*  Copyright PetaChem, LLC 2008-                                      */
/*  No portion may be distributed, modified, used, or divulged without */
/*  express written permission from the copyright holders.             */
/*                                                                     */
/***********************************************************************/
#include "mpi.h"
#include <cstdlib>
#include <math.h>
#define MAX_DATA 200000
#include <cstring>

using namespace std;

int main(int argc, char* argv[])
{

  int UseMPI=1;
  int MPI_SCF_DIE=0;
  char terachem_port_name[1024]; 
  strcpy(terachem_port_name,"terachem_port");

  // DH: added initializations
  int MMAtom;
  int Dx, Dy, Dz;
  int Dxmm, Dymm, Dzmm;
  int totNumAtoms;
  int MMtotNumAtoms;
  int Mode;

  MPI_Comm amber_client;
  MPI_Status mpi_status;
  char port_name[MPI_MAX_PORT_NAME];
  char   buftype[MAX_DATA];
  char   bufchar[100*128];
  char   tmpbuf[1024];
  double bufcoords[MAX_DATA];
  double bufchrgs[MAX_DATA];
  int    bufints[MAX_DATA];
  int mpi_tag = 0;

  // Determine terachem MPI named port name. If already defined through 
  // the CLI, we will respect that.  Alternatively, it may be provided 
  // as the second undashed option argument. Default was set above 
  // (terachem_port)
  if ((argc > 2) && (strcmp(terachem_port_name,"terachem_port")==0)) {
    printf("Assigning ID to terachem_port_name:  %s \n", argv[2]);
    strcpy(terachem_port_name,argv[2]);
  }
  
  int size_mpi, again_mpi;
  int number  = 1;
  MPI_Init(0, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &size_mpi);
  
  // MPI_INFO_NULL are the implementation defaults. 
  MPI_Open_port(MPI_INFO_NULL, port_name);  
  // establishes a port at which the server may be contacted.  
  printf("\nTerachem server available at port_name:   %s  \n",port_name);
  printf("Will publish this port_name to '%s'. \n", terachem_port_name);
  
  // Publish the port name 
  MPI_Publish_name(terachem_port_name, MPI_INFO_NULL, port_name);
  
  printf("Waiting to accept MPI communication from AMBER client. \n");
  MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &amber_client );  // accept communication from AMBER client. 
  printf("MPI communication accepted.  \n");    

  int loop = 1;
  int loop_counter = 0;
// Will go through this loop until MPI gives an exit signal. 
// If !UseMPI, only go through the loop once  
  while (loop) {  
      printf("\nReceiving new QM data.  \n");
      //  Read in num QM atoms, QM atom types, QM atom coordinates
      MPI_Recv( bufints, MAX_DATA, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, amber_client, &mpi_status );
      printf("\nReceived new QM data.  \n");
      if(mpi_status.MPI_TAG == 0) {
	printf("Got a 0 tag from client. Freeing and finalizing.  \n");
	MPI_Comm_free( &amber_client );
	MPI_Unpublish_name(terachem_port_name, MPI_INFO_NULL, port_name);
	MPI_Close_port(port_name);
	MPI_Finalize();
	loop = 0;
        // DH: This exit is potentially dangerous,
        // because statements after while loop are not executed. 
        // This should probably be "break"
	break;
      }
      if(mpi_status.MPI_TAG == 1) {
	printf("Got a 1 tag from client --> Energy mode.  \n");
	Mode = 0; 
      }
      if(mpi_status.MPI_TAG == 2) {
	printf("Got a 2 tag from client --> Gradient mode.  \n");
	Mode = 1;
      }
      
      // atom types - these are stored with two characters each 
      // in the buftype character array
      MPI_Recv( buftype, MAX_DATA, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, 
								amber_client, &mpi_status );
      // Recieve QM coordinates from AMBER
      //printf("Receiving QM coordinates:   \n");
      MPI_Recv( bufcoords, MAX_DATA, MPI_DOUBLE, MPI_ANY_SOURCE, 
								MPI_ANY_TAG, amber_client, &mpi_status );
    // QMMM STARTUP 
    
    // QMMM not yet implemented in ABIN
    bool QMMM = false;
    MMtotNumAtoms = 0;
    MMAtom = 0;
    if(QMMM) {
	printf("Receiving new MM data.  \n");
	// Receive number of MM point charges 
	MPI_Recv( bufints, MAX_DATA, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
		  amber_client, &mpi_status );
	MMtotNumAtoms = bufints[0];
	printf("Number of MM point charges: %d\n", MMtotNumAtoms);
	// Receive new MM charges 
	MPI_Recv( bufchrgs, MAX_DATA, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
		  amber_client, &mpi_status );
	double MMcharge = 0.0;
	// Recieve new coordinates from ABIN
	// printf("\n****** MM coordinates ******\n");
	MPI_Recv( bufcoords, MAX_DATA, MPI_DOUBLE, MPI_ANY_SOURCE, 
		  MPI_ANY_TAG, amber_client, &mpi_status );
    }
    
    // Born not implemented in ABIN
    bool RunBorn = false;
    if(RunBorn) {
      printf("Receiving Born radii.  \n");
      double born_radii[totNumAtoms + MMtotNumAtoms];
      MPI_Recv( bufcoords, MAX_DATA, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, amber_client, &mpi_status );
      for(int atom = 0; atom < totNumAtoms + MMtotNumAtoms; atom ++) {
	born_radii[atom] = bufcoords[atom];  
	printf("Born radii for atom %d:  %lf \n", atom+1, born_radii[atom]); 
      } 
    } 
    
    
   double SCFEnergy; 

   // We need to implement this in ABIN.
   bool MPI_SCF_DIE = false;

    SCFEnergy = 1.0+loop_counter;
    bufcoords[0] = SCFEnergy;
    if(MPI_SCF_DIE) {
        mpi_tag = 1; // send the error message
        printf("SCF did not converge. Setting MPI_TAG = 1. \n");
        MPI_Send( bufcoords, 1, MPI_DOUBLE, 0, mpi_tag, amber_client );
        MPI_SCF_DIE = 0;  // reset the flag for the next loop
        mpi_tag = 0; 
    } else {
        printf("Sending QM energy, QM population charges, and dipoles (QM, MM and total) via MPI. \n");
        printf("QM energy = %.8f \n" , SCFEnergy);
        // Send the energy
        MPI_Send( bufcoords, 1, MPI_DOUBLE, 0, mpi_tag, amber_client );
        // Compute the population charges 
        for(int atom = 0; atom < totNumAtoms; atom ++) {
          bufcoords[atom] = -1-atom; 
        } 
        MPI_Send( bufcoords, totNumAtoms, MPI_DOUBLE, 0, mpi_tag, amber_client );
        // QM dipole moment  
        Dx = -0.01 + loop_counter;
        Dy = -0.02 + loop_counter;
        Dz = -0.03 + loop_counter;
        bufcoords[0] = Dx; 
        bufcoords[1] = Dy; 
        bufcoords[2] = Dz;
        bufcoords[3] = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
        //printf (" DIPOLE: %lf %lf %lf %lf \n", bufcoords[0], bufcoords[1], bufcoords[2], bufcoords[3] ); 
        MPI_Send( bufcoords, 4, MPI_DOUBLE, 0, mpi_tag, amber_client );
        // MM dipole moment 
        Dxmm = -0.01 + loop_counter;
        Dymm = -0.02 + loop_counter;
        Dzmm = -0.03 + loop_counter;
        bufcoords[0] = Dxmm;
        bufcoords[1] = Dymm;
        bufcoords[2] = Dzmm; 
        bufcoords[3] = sqrt(Dxmm*Dxmm + Dymm*Dymm + Dzmm*Dzmm);                  
        //printf (" DIPOLE: %lf %lf %lf %lf \n", bufcoords[0], bufcoords[1], bufcoords[2], bufcoords[3] ); 
        MPI_Send( bufcoords, 4, MPI_DOUBLE, 0, mpi_tag, amber_client );
        // Total dipole moment 
        Dx += Dxmm;
        Dy += Dymm; 
        Dz += Dzmm; 
        bufcoords[0] = Dx;
        bufcoords[1] = Dy;
        bufcoords[2] = Dz; 
        bufcoords[3] = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);                
        //printf (" DIPOLE: %lf %lf %lf %lf \n", bufcoords[0], bufcoords[1], bufcoords[2], bufcoords[3] );  
        MPI_Send( bufcoords, 4, MPI_DOUBLE, 0, mpi_tag, amber_client );
    } 
    if(UseMPI && Mode == 1) {
      printf("Sending gradients via MPI. \n");
      for(int i = 0; i < (totNumAtoms); i ++) {
        bufcoords[3*i]   = 0.001+0.0001*(3*i);
        bufcoords[3*i+1] = 0.001+0.0001*(3*i+1);
        bufcoords[3*i+2] = 0.001+0.0001*(3*i+2);
      }
      for(int i = 0; i < (MMtotNumAtoms); i ++) {
        bufcoords[3*totNumAtoms + 3*i  ] = 0.0;
        bufcoords[3*totNumAtoms + 3*i+1] = 0.0;
        bufcoords[3*totNumAtoms + 3*i+2] = 0.0;
      }
      MPI_Send( bufcoords, 3*totNumAtoms + 3*MMtotNumAtoms, 
      	  MPI_DOUBLE, 0, mpi_tag, amber_client );
    }
      
    loop_counter++; 
    if(!UseMPI)
      loop = 0;
  }  // end while loop for dynamics via MPI
  
  return(0);
}

