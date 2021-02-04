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

    int status = tc.receive();
    // TODO: Call receive individually and validate received data.
    if (status == MPI_TAG_EXIT) {
      break;
    }

    tc.send(loop_counter);
      
    // This is just a precaution, we don't want endless loop!
    loop_counter++;
    if (loop_counter > MAX_LOOP_COUNT)
      break;
  }
  
  return(0);
}
