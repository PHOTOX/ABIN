#include <cstdlib>
#include <math.h>
#include <cstring>
#include "mpi.h"

#define MAX_DATA 200

// When we get an error tag from client, we exit at any time
#define Check_MPI_Recv(stat) \
    if (stat.MPI_TAG == MPI_TAG_ERROR) {\
    throw "Client sent an error tag.";}

using namespace std;

// Buffers for MPI_Recv and MPI_Send calls
char   bufchars[MAX_DATA];
double bufdoubles[MAX_DATA];
int    bufints[MAX_DATA];

static const int MPI_TAG_EXIT = 0;
static const int MPI_TAG_GRADIENT = 2;
static const int MPI_TAG_ERROR = 13;

// This is sent to ABIN for normal MPI_Send
// it is confusing that this is not symmetric, i.e.
// 0 from TC means OK, while 0 from ABIN means exit.
static const int MPI_TAG_OK = 0;

static const int MPI_TAG_SCF_DIED = 1;

// TODO: Separate this to a separate file
// so that it can be re-used in multiple tests.
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

    ~TCServerMock(void) {
      printf("Freeing and finalizing MPI.\n");
      MPI_Comm_free(&abin_client);
      MPI_Close_port(port_name);
      MPI_Finalize();
    }

    // TODO for now we return the client
    //void initializeCommunication() {
    void initializeCommunication() {
      // MPI_INFO_NULL are the implementation defaults. 
      MPI_Open_port(MPI_INFO_NULL, port_name);
      // establishes a port at which the server may be contacted.
      printf("\nTerachem server available at port_name: %s\n", port_name);
      printf("Will publish this port_name to '%s'.\n", terachem_port_name);
      
      // Publish the port name 
      MPI_Publish_name(terachem_port_name, MPI_INFO_NULL, port_name);
      printf("Waiting to accept MPI communication from ABIN client.\n");
      fflush(stdout);

      MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &abin_client);
      printf("MPI communication accepted.\n");

      // It's important to Unpublish the port_name early, otherwise
      // we could get conflicts when other server tried to use the same name.
      MPI_Unpublish_name(terachem_port_name, MPI_INFO_NULL, port_name);
    }

    // This is called only once at the beginning.
    int receiveNumAtoms() {
      printf("Receiving number of atoms...\n");
      fflush(stdout);
      MPI_Recv(bufints, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpi_status);
      Check_MPI_Recv(mpi_status);
      totNumAtoms = bufints[0];
      if (totNumAtoms < 1) {
        printf("ERROR: Invalid number of atoms. Expected positive number, got: %d", totNumAtoms);
        throw "Invalid number of atoms";
      }
      printf("totNumAtoms=%d\n", totNumAtoms);
      return totNumAtoms;
    }

    void receiveAtomTypes() {
      printf("Receiving atom types...\n");
      fflush(stdout);
      MPI_Recv(bufchars, totNumAtoms*2, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpi_status);
      Check_MPI_Recv(mpi_status);
      puts(bufchars);
    }

    // TODO: Validated this!
    void receiveAtomTypesAndScrdir() {
      printf("Receiving atom types and scrdir...\n");
      fflush(stdout);
      MPI_Recv(bufchars, MAX_DATA, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpi_status);
      Check_MPI_Recv(mpi_status);
      // TODO: parse and validate scrdir name.
      // This is a horrible hack in TC
      // so that ABIN can change scratch directories for different beads in PIMD,
      // while retaining the existing Amber interface.
      puts(bufchars);
    }

    // Receive number of QM atoms, QM atom types, QM atom coordinates
    // This is called repeatedly in the MD loop.
    int receive() {
        printf("\nReceiving new QM data.\n");
        fflush(stdout);
        MPI_Recv(bufints, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpi_status);
        Check_MPI_Recv(mpi_status);

        int tag = mpi_status.MPI_TAG;
 
        if (tag == MPI_TAG_EXIT) {
          printf("Got an EXIT tag from client. \n");
          return tag;
        } else if (tag != MPI_TAG_GRADIENT) {
          printf("ERROR: Expected mpi_tag=%d, got %d", MPI_TAG_GRADIENT, tag);
          throw "Invalid MPI TAG received";
        }
        
        if (totNumAtoms != bufints[0]) {
          printf("ERROR: Unexpected number of atoms.\n");
          printf("Expected %d, got %d\n", totNumAtoms, bufints[0]);
          throw "Invalid number of atoms received";
        }

        // TODO: Check that we get same atom types
        // every iteration!
        receiveAtomTypesAndScrdir();
 
        // Receive QM coordinates from ABIN
        // TODO: Separate this to a function.
        printf("Receiving QM coordinates...\n");
        MPI_Recv(bufdoubles, totNumAtoms * 3, MPI_DOUBLE, MPI_ANY_SOURCE, 
            	MPI_ANY_TAG, abin_client, &mpi_status);
        for (int i = 0; i < totNumAtoms * 3; i++) {
          printf("%g ", bufdoubles[i]);
        }
        printf("\n");

        return tag;
    }

    void send(int loop_counter) {
      double SCFEnergy = 1.0 + loop_counter;
      bufdoubles[0] = SCFEnergy;
      int MPI_SCF_DIE = 0;
      if (MPI_SCF_DIE) {
          // TODO: Actually test this scenario!
          printf("SCF did not converge. Setting MPI_TAG_OK = %d.\n", MPI_TAG_SCF_DIED);
          MPI_Send(bufdoubles, 1, MPI_DOUBLE, 0, MPI_TAG_SCF_DIED, abin_client);
          MPI_SCF_DIE = 0; // reset the flag for the next loop
      }
     
      printf("Sending QM energy, QM population charges, and dipoles (QM, MM and total) via MPI.\n");
      printf("QM energy = %.8f \n" , SCFEnergy);
      // Send the energy
      MPI_Send(bufdoubles, 1, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client );
      // Compute the population charges 
      for(int atom = 0; atom < totNumAtoms; atom++) {
        bufdoubles[atom] = -1 - atom; 
      }
      MPI_Send(bufdoubles, totNumAtoms, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
      // QM dipole moment
      double Dx = -0.01 + loop_counter;
      double Dy = -0.02 + loop_counter;
      double Dz = -0.03 + loop_counter;
      double DTotal = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
      bufdoubles[0] = Dx;
      bufdoubles[1] = Dy;
      bufdoubles[2] = Dz;
      bufdoubles[3] = DTotal;
      printf("QM DIPOLE: %lf %lf %lf %lf\n", Dx, Dy, Dz, DTotal);
      MPI_Send(bufdoubles, 4, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
     
      // NOTE: In the real TC interface, gradients are sent
      // conditionally only if they are requested.
      // But we always request them in ABIN (see tc.receive).
      // TODO: Maybe we could link this to the water force field
      // in waterpotentials/ and get real forces.
      printf("Sending gradients via MPI. \n");
      for(int i = 0; i < (totNumAtoms); i ++) {
        bufdoubles[3*i]   = 0.001+0.0001*(3*i);
        bufdoubles[3*i+1] = 0.001+0.0001*(3*i+1);
        bufdoubles[3*i+2] = 0.001+0.0001*(3*i+2);
      }
      MPI_Send(bufdoubles, 3*totNumAtoms, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
    }


  private:
    char terachem_port_name[1024];
    char port_name[MPI_MAX_PORT_NAME];
    MPI_Comm abin_client;
    MPI_Status mpi_status;
    int totNumAtoms;
};

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
