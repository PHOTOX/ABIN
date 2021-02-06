#include "tc_mpi_api.h"

using namespace std;

TCServerMock::TCServerMock(char *portname) {
  strcpy(terachemPortName, portname);
  // Initialize MPI in the constructor
  MPI_Init(0, NULL);

  // Check that we're only running one MPI process
  // i.e. that we've been invoked by `mpirun -n 1`
  int commSize;
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  printf("MPI_Comm_size = %d\n", commSize);
  if (commSize != 1) {
    printf("ERROR: Comm_size != 1!\n");
    printf("Please execute this program with 'mpirun -n 1'\n");
    throw "Incorrect mpirun invocation";
  }
}

TCServerMock::~TCServerMock(void) {
  printf("Freeing and finalizing MPI.\n");
  MPI_Comm_free(&abin_client);
  MPI_Close_port(mpiPortName);
  MPI_Finalize();
}

void TCServerMock::checkRecvCount(MPI_Status *mpiStatus,
                                  MPI_Datatype datatype,
                                  int expected_count) {
  int recvCount;
  MPI_Get_count(mpiStatus, datatype, &recvCount);
  if (recvCount != expected_count) {
    printf("Unexpected received count\n");
    printf("Expected %d, got %d\n", expected_count, recvCount);
    throw "Unexpected received count";
  }
}

void TCServerMock::initializeCommunication() {
  // MPI_INFO_NULL are the implementation defaults. 
  MPI_Open_port(MPI_INFO_NULL, mpiPortName);
  // establishes a port at which the server may be contacted.
  printf("Terachem server available at port name: %s\n", mpiPortName);
  printf("Will publish this port under name '%s' to hydra_nameserver.\n", terachemPortName);
  
  // Publish the port name 
  MPI_Publish_name(terachemPortName, MPI_INFO_NULL, mpiPortName);
  printf("Waiting to accept MPI communication from ABIN client.\n");
  fflush(stdout);

  MPI_Comm_accept(mpiPortName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &abin_client);
  printf("MPI communication accepted.\n");

  // It's important to Unpublish the port_name early, otherwise
  // we could get conflicts when other server tried to use the same name.
  MPI_Unpublish_name(terachemPortName, MPI_INFO_NULL, mpiPortName);
}

// This is called only once at the beginning.
int TCServerMock::receiveNumAtoms() {
  printf("Receiving number of atoms...\n");
  fflush(stdout);
  int recvCount = 1;
  MPI_Recv(bufints, recvCount, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_INT, recvCount);
  totNumAtoms = bufints[0];
  if (totNumAtoms < 1) {
    printf("ERROR: Invalid number of atoms. Expected positive number, got: %d", totNumAtoms);
    throw "Invalid number of atoms";
  }
  printf("totNumAtoms=%d\n", totNumAtoms);
  return totNumAtoms;
}

// TODO: Each receive function should actually return 
// the data it received so that they can be validated.
void TCServerMock::receiveAtomTypes() {
  printf("Receiving atom types...\n");
  fflush(stdout);
  int recvCount = totNumAtoms * 2;
  MPI_Recv(bufchars, recvCount, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_CHAR, recvCount);
  puts(bufchars);
}

void TCServerMock::receiveAtomTypesAndScrdir() {
  // TODO: Check that we get same atom types every iteration!
  printf("Receiving atom types and scrdir...\n");
  fflush(stdout);
  MPI_Recv(bufchars, MAX_DATA, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  // TODO: parse and validate scrdir name.
  // This is a horrible hack in TC
  // so that ABIN can change scratch directories for different beads in PIMD,
  // while retaining the existing Amber interface.
  puts(bufchars);
}

void TCServerMock::receiveCoordinates() {
  // Receive QM coordinates from ABIN
  // TODO: Separate this to a function.
  printf("Receiving QM coordinates...\n");
  int recvCount = totNumAtoms * 3;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE,
          MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  for (int i = 0; i < totNumAtoms * 3; i++) {
    printf("%g ", bufdoubles[i]);
  }
  printf("\n");
}

// Receive number of QM atoms, QM atom types, QM atom coordinates
// This is called repeatedly in the MD loop.
int TCServerMock::receiveBeginLoop() {
  printf("\nReceiving new QM data.\n");
  fflush(stdout);
  int recvCount = 1;
  MPI_Recv(bufints, recvCount, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  // When ABIN sends it's exit tag, we want to allow it
  // to send zero data.
  //checkRecvCount(&mpiStatus, MPI_INT, recvCount);
  int tag = mpiStatus.MPI_TAG;

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
  return tag;
}

// This is what we expect from ABIN every MD iteration
// for classical MD and PIMD.
// For Surface Hopping see receive_sh()
int TCServerMock::receive() {
  int tag = receiveBeginLoop();
  if (tag == MPI_TAG_EXIT) {
    return tag;
  }
  receiveAtomTypesAndScrdir();
  receiveCoordinates();
  return tag;
}

void TCServerMock::sendSCFEnergy(double energy, int MPI_SCF_DIE) {
  bufdoubles[0] = energy;
  if (MPI_SCF_DIE) {
    // TODO: Actually test this scenario!
    printf("SCF did not converge. Setting MPI_TAG_OK = %d.\n", MPI_TAG_SCF_DIED);
    MPI_Send(bufdoubles, 1, MPI_DOUBLE, 0, MPI_TAG_SCF_DIED, abin_client);
    MPI_SCF_DIE = 0; // reset the flag for the next loop
  }
 
  printf("QM energy = %.8f\n" , energy);
  MPI_Send(bufdoubles, 1, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

void TCServerMock::sendQMCharges() {
  // TODO: We should move this computation elsewhere
  // and accept charges as inputs here.
  //
  // Compute fake population charges
  for(int atom = 0; atom < totNumAtoms; atom++) {
    bufdoubles[atom] = -1 - atom; 
  }
  MPI_Send(bufdoubles, totNumAtoms, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

void TCServerMock::sendQMDipoleMoments() {
  // QM dipole moment
  double Dx = -0.01;
  double Dy = -0.02;
  double Dz = -0.03;
  double DTotal = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
  bufdoubles[0] = Dx;
  bufdoubles[1] = Dy;
  bufdoubles[2] = Dz;
  bufdoubles[3] = DTotal;
  printf("QM DIPOLE: %lf %lf %lf %lf\n", Dx, Dy, Dz, DTotal);
  MPI_Send(bufdoubles, 4, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

void TCServerMock::sendQMGradients() {
  printf("Sending gradients via MPI. \n");
  for(int i = 0; i < (totNumAtoms); i ++) {
    bufdoubles[3*i]   = 0.001+0.0001*(3*i);
    bufdoubles[3*i+1] = 0.001+0.0001*(3*i+1);
    bufdoubles[3*i+2] = 0.001+0.0001*(3*i+2);
    printf("%lf %lf %lf\n", bufdoubles[3*i], bufdoubles[3*i+1], bufdoubles[3*i+2]);
  }
  MPI_Send(bufdoubles, 3*totNumAtoms, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

void TCServerMock::send(int loop_counter) {
  printf("Sending QM energy, QM population charges, QM dipoles and QN gradients via MPI.\n");

  // TODO: Maybe we could link this to the water force field
  // in waterpotentials/ and get real energies and forces?
  double SCFEnergy = 1.0 + loop_counter;
  int MPI_SCF_DIE = 0;
  sendSCFEnergy(SCFEnergy, MPI_SCF_DIE);

  sendQMCharges();

  sendQMDipoleMoments();

  // NOTE: In the real TC interface, gradients are sent
  // conditionally only if they are requested.
  // But we always request them in ABIN (see tc.receive).
  sendQMGradients();
}
