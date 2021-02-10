#include "tc_mpi_api.h"

using namespace std;

TCServerMock::TCServerMock(char *serverName) {

  tcServerName = NULL;
  gradients = coordinates = NULL;
  if (serverName) {
    tcServerName = new char[1024];
    strcpy(tcServerName, serverName);
  }

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
    throw std::runtime_error("Incorrect mpirun invocation");
  }
}

TCServerMock::~TCServerMock(void) {
  printf("Freeing and finalizing MPI.\n");
  MPI_Comm_free(&abin_client);
  MPI_Close_port(mpiPortName);
  MPI_Finalize();
  if (gradients) {
    delete[] gradients;
  }
  if (coordinates) {
    delete[] coordinates;
  }
}

void TCServerMock::checkRecvCount(MPI_Status *mpiStatus,
                                  MPI_Datatype datatype,
                                  int expected_count) {
  int recvCount;
  MPI_Get_count(mpiStatus, datatype, &recvCount);
  if (recvCount != expected_count) {
    printf("Unexpected received count\n");
    printf("Expected %d, got %d\n", expected_count, recvCount);
    throw std::runtime_error("Unexpected received count");
  }
}

// When we get an error tag from client we exit immediately.
void TCServerMock::checkRecvTag(MPI_Status &mpiStatus) {
  if (mpiStatus.MPI_TAG == MPI_TAG_ERROR) {
    throw std::runtime_error("Client sent an error tag.");
  }
}

void TCServerMock::initializeCommunication() {
  // MPI_INFO_NULL are the implementation defaults. 
  MPI_Open_port(MPI_INFO_NULL, mpiPortName);
  // Establishes a port at which the server may be contacted.
  printf("Fake TeraChem server available at port name: %s\n", mpiPortName);
  
  // Publish the port name, but only if tcServerName was passed to constructor.
  if (tcServerName) {
    MPI_Publish_name(tcServerName, MPI_INFO_NULL, mpiPortName);
    printf("Port published under server name '%s'\n", tcServerName);
  } else {
    printf("Server name not specified, nothing to publish\n");
    printf("Pass the port name manually.\n");
  }

  printf("Waiting to accept MPI communication from ABIN client.\n");
  fflush(stdout);

  MPI_Comm_accept(mpiPortName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &abin_client);
  printf("MPI communication accepted.\n");

  // It's important to unpublish the port_name early, otherwise
  // we could get conflicts when other server tried to use the same name.
  // WARNING: If more then one tc_server is running,
  // concurrent calls MPI_Unpublish_name are crashing the hydra_nameserver.
  // https://github.com/pmodels/mpich/issues/5058
  if (tcServerName) {
    MPI_Unpublish_name(tcServerName, MPI_INFO_NULL, mpiPortName);
  }

  delete[] tcServerName;
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
    throw std::runtime_error("Invalid number of atoms");
  }
  printf("totNumAtoms=%d\n", totNumAtoms);

  // Allocate buffers for gradients and coordinates
  gradients = new double[totNumAtoms * 3];
  coordinates = new double[totNumAtoms * 3];
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
  // TODO: Check that we get the same atom types every iteration!
  printf("Receiving atom types and scrdir...\n");
  fflush(stdout);
  MPI_Recv(bufchars, MAX_DATA, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  // TODO: parse and validate scrdir name.
  // This is a horrible hack in TC
  // ABIN misuses the atom type array to set scratch directories
  // (useful e.g. for different beads in PIMD)
  // This was done this way to preserve the existing Amber interface.
  puts(bufchars);
}

void TCServerMock::receiveCoordinates() {
  // Receive QM coordinates from ABIN
  printf("Receiving QM coordinates...\n");
  int recvCount = totNumAtoms * 3;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE,
          MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  for (int i = 0; i < totNumAtoms * 3; i++) {
    coordinates[i] = bufdoubles[i];
    printf("%g ", bufdoubles[i]);
  }
  printf("\n");
}

// Receive number of QM atoms, QM atom types, QM atom coordinates.
// This is called repeatedly in the MD loop.
int TCServerMock::receiveBeginLoop() {
  printf("\nReceiving new QM data.\n");
  fflush(stdout);
  int recvCount = 1;
  MPI_Recv(bufints, recvCount, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  // When ABIN sends it's exit tag,
  // we want to allow it to send zero data.
  //checkRecvCount(&mpiStatus, MPI_INT, recvCount);
  int tag = mpiStatus.MPI_TAG;

  if (tag == MPI_TAG_EXIT) {
    printf("Got an EXIT tag from client. \n");
    return tag;
  } else if (tag != MPI_TAG_GRADIENT) {
    printf("ERROR: Expected mpi_tag=%d, got %d", MPI_TAG_GRADIENT, tag);
    throw std::runtime_error("Invalid MPI TAG received");
  }
  
  if (totNumAtoms != bufints[0]) {
    printf("ERROR: Unexpected number of atoms.\n");
    printf("Expected %d, got %d\n", totNumAtoms, bufints[0]);
    throw std::runtime_error("Invalid number of atoms received");
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
    printf("SCF did not converge. Setting MPI_TAG = %d.\n", MPI_TAG_SCF_DIED);
    MPI_Send(bufdoubles, 1, MPI_DOUBLE, 0, MPI_TAG_SCF_DIED, abin_client);
    MPI_SCF_DIE = 0; // reset the flag for the next loop
  }
 
  printf("QM energy = %.8f\n" , energy);
  MPI_Send(bufdoubles, 1, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

// Compute fake Mulliken charges.
void TCServerMock::computeFakeQMCharges(double *charges) {
  for (int atom = 0; atom < totNumAtoms; atom++) {
    charges[atom] = -1 - atom; 
  }
}

void TCServerMock::sendQMCharges() {
  // TODO: We could move this computation elsewhere
  // and accept charges as inputs here.
  computeFakeQMCharges(bufdoubles);
  MPI_Send(bufdoubles, totNumAtoms, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

// Compute fake dipole moments.
void TCServerMock::computeFakeQMDipoleMoment(double &Dx, double &Dy, double &Dz, double &DTotal) {
  Dx = -0.01;
  Dy = -0.02;
  Dz = -0.03;
  DTotal = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
  printf("QM DIPOLE: %lf %lf %lf %lf\n", Dx, Dy, Dz, DTotal);
}

void TCServerMock::sendQMDipoleMoment() {
  computeFakeQMDipoleMoment(bufdoubles[0], bufdoubles[1], bufdoubles[2], bufdoubles[3]);
  MPI_Send(bufdoubles, 4, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

void TCServerMock::sendQMGradients() {
  printf("Sending gradients via MPI. \n");
  for(int i = 0; i < totNumAtoms; i++) {
    bufdoubles[3*i]   = gradients[3*i];
    bufdoubles[3*i+1] = gradients[3*i+1];
    bufdoubles[3*i+2] = gradients[3*i+2];
    printf("%lf %lf %lf\n", bufdoubles[3*i], bufdoubles[3*i+1], bufdoubles[3*i+2]);
  }
  MPI_Send(bufdoubles, 3 * totNumAtoms, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
}

void TCServerMock::send() {
  printf("Sending QM energy, QM population charges, QM dipoles and QN gradients via MPI.\n");

  double energy = getWaterGradients();
  int MPI_SCF_DIE = 0;

  sendSCFEnergy(energy, MPI_SCF_DIE);

  sendQMCharges();

  sendQMDipoleMoment();

  // NOTE: In the real TC interface, gradients are sent
  // conditionally only if they are requested.
  // But we always request them in ABIN (see tc.receive).
  sendQMGradients();
}

// Gradients stored internally for now,
// returns energy in atomic units.
double TCServerMock::getWaterGradients() {
  // conversion constants
  const double ANG = 1.889726132873;
  const double AUTOKCAL = 627.50946943;
  const double FORCE_FAC = 1 / ANG / AUTOKCAL;

  double *converted_coords = new double[3 * totNumAtoms];

  if (totNumAtoms % 3 != 0) {
    throw std::runtime_error("Number of atoms not divisible by 3");
  }
  int nwater = totNumAtoms / 3;

  double E = 0.0;
  for (int i = 0; i < totNumAtoms * 3; i++) {
    converted_coords[i] = coordinates[i] / ANG;
    gradients[i] = 0;
  }

  h2o::qtip4pf water_potential;
  E = water_potential(nwater, coordinates, gradients);

  // Convert energy and gradients to atomic units
  E /= AUTOKCAL;
  for (int i = 0; i < totNumAtoms * 3; i++) {
     gradients[i] *= FORCE_FAC;
  }

  delete[] converted_coords;

  return E;
}
