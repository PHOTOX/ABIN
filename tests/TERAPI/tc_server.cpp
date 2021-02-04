#include <cstdlib>
#include <cstring>

#include "../tc_mpi_api.h"

using namespace std;

int main(int argc, char* argv[])
{
  char terachem_port_name[1024]; 

  strcpy(terachem_port_name, "terachem_port.1");
  if (argc > 1) {
    strcpy(terachem_port_name,argv[1]);
  }

  TCServerMock tc = TCServerMock(terachem_port_name);

  tc.initializeCommunication();

  tc.receiveNumAtoms();
  tc.receiveAtomTypes();

  int loop_counter = 0;
  int MAX_LOOP_COUNT = 100;
  // Will go through this loop until MPI client gives an exit signal.
  while (true) {

    int status = tc.receive_begin_loop();
    // TODO: Call receive individually and validate received data.
    if (status == MPI_TAG_EXIT) {
      break;
    }

    // TODO: Validate received data?
    // Should be hardcode what we expect to receive?
    tc.receiveAtomTypesAndScrdir();
    tc.receiveCoordinates();

    //tc.send(loop_counter);

    // TODO: Maybe we could link this to the water force field
    // in waterpotentials/ and get real energies and forces?
    double SCFEnergy = 1.0 + loop_counter;
    int MPI_SCF_DIE = 0;
    tc.send_scf_energy(SCFEnergy, MPI_SCF_DIE);
 
    tc.send_qm_charges();
 
    tc.send_qm_dipole_moments();

    // NOTE: In the real TC interface, gradients are sent
    // conditionally only if they are requested.
    // But we always request them in ABIN (see tc.receive).
    tc.send_qm_gradients();
      
    // This is just a precaution, we don't want endless loop!
    loop_counter++;
    if (loop_counter > MAX_LOOP_COUNT)
      break;
  }
  
  return(0);
}
