#include <cstdlib>
#include <cstring>

#include "../tc_mpi_api.h"

using namespace std;

int main(int argc, char* argv[])
{
  char *serverName = NULL;

  // Due to a bug in hydra_nameserver, it crashes
  // when multiple TC servers call `MPI_Unpublish_name()`
  // Hence, we want to allow invoking without this parameter,
  // in which case TC server will just print the port to stdin,
  // where it could be grepped and passed via file to ABIN,
  // and it will never call MPI_Publish_name/MPI_Unpublish_name
  // NOTE: This behaviour is different from real TC,
  // which has default serverName and will always try to publish it.
  if (argc > 2) {
    printf("Only one cmdline argument supported, <serverName>, but you provided more!");
    throw std::runtime_error("Incorrect invocation");
  }

  if (argc == 2) {
    serverName = new char[1024];
    strcpy(serverName, argv[1]);
  }

  TCServerMock tc = TCServerMock(serverName);

  if (serverName) {
    delete[] serverName;
  }

  tc.initializeCommunication();

  tc.receiveNumAtoms();
  tc.receiveAtomTypes();

  int status = tc.receive();
  if (status == MPI_TAG_EXIT) {
    throw std::runtime_error("unexpected exit tag");
  }

  // At this point ABIN expects SCF energy. Let's send
  // zero doubles instead of one to throw it off!
  printf("Sending invalid (zero bytes) data, muhehehe!\n");
  MPI_Comm *abinComm = tc.getABINCommunicator();
  double energies[1] = {1};
  MPI_Send(energies, 0, MPI_DOUBLE, 0, 0, *abinComm);

  // Calling receive since we're expecting MPI_ERROR_TAG from ABIN.
  tc.receive();
  return(0);
}
