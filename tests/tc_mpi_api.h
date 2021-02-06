#include <cstdlib>
#include <math.h>
#include <cstring>
#include "mpi.h"

#define MAX_DATA 200

// When we get an error tag from client, we exit at any time
#define checkRecvTag(mpi_status) \
    if (mpi_status.MPI_TAG == MPI_TAG_ERROR) {\
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
    TCServerMock(char *);
    ~TCServerMock(void);

    void initializeCommunication();

    int receiveNumAtoms();
    void receiveAtomTypes();
    void receiveAtomTypesAndScrdir();
    void receiveCoordinates();

    int receiveBeginLoop();
    // This combines all expected receive
    // methods in their order. But the tc_server
    // implementation can also call them individually
    // for better granurality.
    int receive();

    void send(int loop_counter);
    void sendSCFEnergy(double, int);
    void sendQMCharges();
    void sendQMDipoleMoments();
    void sendQMGradients();

  private:
    // This one will be published via MPI_Publish
    // ABIN will be looking for this one via hydra_nameserver
    char terachemPortName[1024];
    // This is the name of the actual port from MPI_Open_port()
    char mpiPortName[MPI_MAX_PORT_NAME];
    // Buffers for MPI_Recv and MPI_Send calls
    char   bufchars[MAX_DATA];
    double bufdoubles[MAX_DATA];
    int    bufints[MAX_DATA];

    MPI_Comm abin_client;
    MPI_Status mpiStatus;

    int totNumAtoms;
    char *atomTypes[MAX_DATA];

    void checkRecvCount(MPI_Status*, MPI_Datatype, int);
};
