#include <cstdlib>
#include <cstring>

#include "../tc_mpi_api.h"

using namespace std;

int main(int argc, char* argv[])
{
  char server_name[1024];

  if (argc != 2) {
    printf("I need exactly one cmdline argument <server_name>");
    throw std::runtime_error("Incorrect invocation");
  }

  strcpy(server_name, argv[1]);

  TCServerMock tc = TCServerMock(server_name);

  tc.initializeCommunication();

  tc.receiveNumAtoms();
  tc.receiveAtomTypes();

  int loop_counter = 0;
  int MAX_LOOP_COUNT = 100;
  // Will go through this loop until MPI client gives an exit signal.
  while (true) {

    int status = tc.receiveBeginLoop();
    if (status == MPI_TAG_EXIT) {
      break;
    }

    tc.receiveAtomTypesAndScrdir();
    tc.receiveCoordinates();

    // Energies and gradients from qTIP4PF potential
    double energy = tc.getWaterGradients();

    int MPI_SCF_DIE = 0;
    tc.sendSCFEnergy(energy, MPI_SCF_DIE);
 
    tc.sendQMCharges();
 
    tc.sendQMDipoleMoments();

    // NOTE: In the real TC interface, gradients are sent
    // conditionally only if they are requested.
    // But we always request them in ABIN (see tc.receive).
    tc.sendQMGradients();
      
    // This is just a precaution, we don't want endless loop!
    loop_counter++;
    if (loop_counter > MAX_LOOP_COUNT) {
      printf("Maximum number of steps exceeded!\n");
      return(1);
    }
  }
  
  return(0);
}
