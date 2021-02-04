#include <cstdlib>
#include <math.h>
#include <cstring>
#include "mpi.h"

#define MAX_DATA 200

// When we get an error tag from client, we exit at any time
#define Check_MPI_Recv(stat) \
    if (stat.MPI_TAG == MPI_TAG_ERROR) {\
    throw "Client sent an error tag.";}

#define MPI_TAG_EXIT 0
#define MPI_TAG_GRADIENT 2
#define MPI_TAG_ERROR 13

// This is sent to ABIN for normal MPI_Send
// it is confusing that this is not symmetric, i.e.
// 0 from TC means OK, while 0 from ABIN means exit.
#define MPI_TAG_OK 0

#define MPI_TAG_SCF_DIED 1

class TCServerMock {
  public:
    TCServerMock(char *portname) {
      strcpy(terachem_port_name, portname);
      // Initialize MPI in the constructor
      // TODO: Check for errors.
      MPI_Init(0, NULL);

      // Check that we're only running one MPI process
      // i.e. that we've been invoked by `mpirun -n 1`
      int mpi_comm_size;
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
      printf("MPI_Comm_size = %d\n", mpi_comm_size);
      if (mpi_comm_size != 1) {
        printf("ERROR: Comm_size != 1!\n");
        printf("Please execute this program with mpirun -n 1\n");
      throw "Incorrect mpirun invocation";
      }
    }

    ~TCServerMock(void);

    void initializeCommunication();

    int receiveNumAtoms();
    void receiveAtomTypes();
    // TODO: Validated this!
    void receiveAtomTypesAndScrdir();
    void receiveCoordinates();

    // Receive number of QM atoms, QM atom types, QM atom coordinates
    // This is called repeatedly in the MD loop.
    int receive();
    int receive_sh();

    void send(int loop_counter);
    void send_sh(int loop_counter);

  private:
    char terachem_port_name[1024];
    char port_name[MPI_MAX_PORT_NAME];
    // Buffers for MPI_Recv and MPI_Send calls
    char   bufchars[MAX_DATA];
    double bufdoubles[MAX_DATA];
    int    bufints[MAX_DATA];

    MPI_Comm abin_client;
    MPI_Status mpi_status;
    int totNumAtoms;
};
