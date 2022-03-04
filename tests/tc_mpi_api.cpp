#include<unistd.h>
#include<fstream>
#include "tc_mpi_api.h"

using namespace std;

namespace TCTensor {
  double *double1D(unsigned long long length) {
    if (length == 0)
      return NULL;

    double *A = (double*) malloc(sizeof(double) * length);
    if (A == NULL) {
      throw std::runtime_error("allocation failed");
    }

    std::memset((void*) A, 0, sizeof(double) * length);
    return A;
  }

  void freeDouble1D(double *&A) {
    if (A == NULL)
      return;
    free(A);
    A = NULL;
  }
}

TCServerMock::TCServerMock(char *serverName) {

  tcServerName = NULL;
  gradients = coordinates = NULL;
  atomTypes = NULL;
  if (serverName) {
    tcServerName = new char[strlen(serverName)+1];
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

  // This needs to be set manually in tc_server.cpp
  FMSNumStates = -1;
}


TCServerMock::~TCServerMock(void) {
  printf("Freeing and finalizing MPI.\n");
  MPI_Comm_free(&abin_client);
  if (mpiPortName) {
    MPI_Close_port(mpiPortName);
  }
  MPI_Finalize();
  if (atomTypes) {
    delete[] atomTypes;
  }
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


void TCServerMock::printMPIError(int error_code) {
    int *resultLen = NULL;
    int new_error_code = MPI_Error_string(error_code, bufchars, resultLen);
    if (new_error_code == MPI_SUCCESS) {
      printf("%s\n", bufchars);
    }
}


MPI_Comm* TCServerMock::getABINCommunicator() {
  return &abin_client;
}


// Publish server name, but only if tcServerName was passed to constructor.
void TCServerMock::publishServerName(char *serverName, char *portName) {
  if (!serverName) {
    printf("Server name not specified, nothing to publish\n");
    printf("Pass the port name manually.\n");
    return;
  }
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  // In seconds
  double sleepTime = 0.2;
  double maxTime = 50;
  double elapsedTime = 0;

  // Retry if hydra_nameserver isn't ready.
  int ierr = MPI_Publish_name(serverName, MPI_INFO_NULL, portName);
  while (ierr != MPI_SUCCESS) {
    if (elapsedTime > maxTime) {
      printMPIError(ierr);
      throw std::runtime_error("could not publish server name");
    }
    sleep(sleepTime);
    elapsedTime += sleepTime;
    ierr = MPI_Publish_name(serverName, MPI_INFO_NULL, portName);
  }

  printf("Port published under server name '%s'\n", tcServerName);

  ierr = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  if (ierr != MPI_SUCCESS) {
    printMPIError(ierr);
    throw std::runtime_error("could not set error handler");
  }
}


void TCServerMock::initializeCommunication() {
  // Establishes a port at which the server may be contacted.
  // MPI_INFO_NULL are the implementation defaults. 
  MPI_Open_port(MPI_INFO_NULL, mpiPortName);

  printf("Fake TeraChem server available at port name: %s\n", mpiPortName);
  fflush(stdout);

  publishServerName(tcServerName, mpiPortName);

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
    delete[] tcServerName;
  }
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


char* TCServerMock::receiveAtomTypes() {
  printf("Receiving atom types...\n");
  fflush(stdout);
  int recvCount = totNumAtoms * 2;
  MPI_Recv(bufchars, recvCount, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_CHAR, recvCount);

  atomTypes = new char[totNumAtoms * 2 + 1];
  strncpy(atomTypes, bufchars, totNumAtoms * 2);
  atomTypes[totNumAtoms * 2] = '\0';
  puts(atomTypes);
  return atomTypes;
}


void TCServerMock::receiveAtomTypesAndScrdir() {
  printf("Receiving atom types and scrdir...\n");
  MPI_Recv(bufchars, MAX_DATA, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  puts(bufchars);
  char *atoms = new char[2*totNumAtoms + 1];
  strncpy(atoms, bufchars, totNumAtoms * 2);
  atoms[totNumAtoms * 2] = '\0';
  if (strcmp(atomTypes, atoms) != 0) {
    printf("ERROR: expected '%s', got '%s'\n", atomTypes, atoms);
    throw std::runtime_error("invalid atom types");
  }
  delete[] atoms;

  char scrdir[1024];
  parseScrDir(bufchars, scrdir);
  validateScrDir(scrdir);
  printf("Using scratch directory %s\n", scrdir);
}

void TCServerMock::validateScrDir(char *scrdir) {
  unsigned int beadIdx;
  int ret = sscanf(scrdir, "scrdir%4u", &beadIdx);
  if (ret == 0 || ret == EOF) {
    throw std::runtime_error("invalid scrdir");
  }
  printf("Bead index: %u\n", beadIdx);
  // To make sure all scrdirs were passed correctly,
  // we create empty files with their names so
  // that they can be compared with the reference.
  std::ofstream output(scrdir);
}

void TCServerMock::parseScrDir(char *buffer, char *scrdir) {
  // DH VARIABLE SCRATCH DIRECTORY, useful for Path Integral MD
  // Done a bit awkwardly to preserve AMBER interface
  // The format in bufchars needs to be: AtomTypes++scrdir++
  // i.e. "C H H H H ++scrdir01++
  // Code copied directly from terachem/amber.cpp
  char delim = '+';
  if (buffer[totNumAtoms*2] == delim && buffer[totNumAtoms*2+1] == delim) {
    int i = totNumAtoms * 2 + 2;
    int j = 0;
    while (bufchars[i] != delim) {
      scrdir[j] = bufchars[i];
      i++;j++;
    }
    scrdir[j] = '\0';
  }
}

void TCServerMock::receiveCoordinates() {
  // Receive QM coordinates from ABIN
  printf("Receiving QM coordinates in angstroms\n");
  const double ANG = 1.889726132873;
  int recvCount = totNumAtoms * 3;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE,
          MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  for (int i = 0; i < totNumAtoms * 3; i++) {
    // Conversion to Bohrs
    coordinates[i] = bufdoubles[i] * ANG;
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

  if (totNumAtoms % 3 != 0) {
    throw std::runtime_error("Number of atoms not divisible by 3");
  }

  int nwater = totNumAtoms / 3;
  double E = 0.0;

  // The water potential expects coordinates in angstromgs
  double *converted_coords = new double[3 * totNumAtoms];
  for (int i = 0; i < totNumAtoms * 3; i++) {
    converted_coords[i] = coordinates[i] / ANG;
    gradients[i] = 0;
  }

  h2o::qtip4pf water_potential;
  E = water_potential(nwater, converted_coords, gradients);

  // Convert energy and gradients to atomic units
  E /= AUTOKCAL;
  for (int i = 0; i < totNumAtoms * 3; i++) {
     gradients[i] *= FORCE_FAC;
  }

  delete[] converted_coords;

  return E;
}

int TCServerMock::fmsinit_receive()
{

  printf("Receiving first set of data from ABIN\n");
  int recvCount = 3;
  MPI_Recv(bufints, recvCount, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_INT, recvCount);

  int tag = mpiStatus.MPI_TAG;
  if(tag == MPI_TAG_EXIT) 
  {
    printf("Got an EXIT tag from client.\n");
    return tag;
  }

  int FMSInit = bufints[0];
  FMS_.NAtoms = bufints[1];
  FMS_.MMAtoms = bufints[2];
  printf("FMSInit flag=%d, NAtoms=%d, MMAtoms=%d\n", FMSInit, FMS_.NAtoms, FMS_.MMAtoms);

  if (FMSInit != 1) {
    throw std::runtime_error("Should get FMS init flag here");
  }

  if (mpiStatus.MPI_TAG != 2) {
     throw std::runtime_error("Only FMS mode is allowed in FMS mode");
  }

  allocate_fms_data();

  // Set up basis set on first call from ABIN
  printf("Receiving atom types from ABIN\n");
  receiveAtomTypes();
  // atom types - these are stored with two characters each in the bufchars character array
  /*
  MPI_Recv(bufchars, MAX_DATA, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  if (strlen(bufchars) < (NAtoms + MMAtoms) * 2) {
     throw std::runtime_error("did not receive atom types");
  }
  */

  // Receive first QM coordinates from ABIN
  // We don't actually need to do anythin with them here,
  // in TC this is just to get initial WF, but nothing is passed back to ABIN, see fmsinit_send()
  printf("Receiving coordinates from ABIN\n");
  recvCount = (FMS_.NAtoms + FMS_.MMAtoms) * 3;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  for (int i = 0; i < recvCount; i++) {
    coordinates[i] = bufdoubles[i];
    printf("%g ", bufdoubles[i]);
  }
  return tag;
}

void TCServerMock::fmsinit_send() {
  // This this gets sent after the first receive from ABIN`
  // TODO: Figure out where to hold all this data
  bufints[0] = Data_.nci;
  bufints[1] = Data_.nbf;
  bufints[2] = Data_.nblob;
  MPI_Send(bufints, 3, MPI_INT, 0, MPI_TAG_OK, abin_client);
}

int TCServerMock::fms_receive() {
  FMS_.FMSInit = 1;
  FMS_.DoCoup = 1;

  printf("\n\nReceiving new FMS data...  \n");

  int recvCount = 12;
  MPI_Recv(bufints, recvCount, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  int tag = mpiStatus.MPI_TAG;
  if (tag == MPI_TAG_EXIT)  {
    printf("Got an EXIT tag from client.\n");
    return tag;
  }
  checkRecvCount(&mpiStatus, MPI_INT, recvCount);

  FMS_.FMSInit = bufints[0];
  FMS_.NAtoms = bufints[1];
  FMS_.DoCoup = bufints[2];
  FMS_.TrajID = bufints[3];
  FMS_.Cent1 = bufints[4];
  FMS_.Cent2 = bufints[5];
  FMS_.StateID = bufints[6];
  FMS_.OldWfn = bufints[7];
  FMS_.iState = bufints[8];
  FMS_.jState = bufints[9];

  FMS_.FirstCall = bufints[10];
  FMS_.FMSRestart = bufints[11];

  printf("FMSInit   :  %d\n", bufints[0]);
  printf("NAtoms    :  %d\n", bufints[1]);
  printf("DoCoup    :  %d\n", bufints[2]);
  printf("TrajID    :  %d\n", bufints[3]);
  printf("Cent1     :  %d\n", bufints[4]);
  printf("Cent2     :  %d\n", bufints[5]);
  printf("StateID   :  %d\n", bufints[6]);
  printf("OldWfn    :  %d\n", bufints[7]);
  printf("iState    :  %d\n", bufints[8]);
  printf("jState    :  %d\n", bufints[9]);
  printf("FirstCall :  %d\n", bufints[10]);
  printf("FMSRestart:  %d\n", bufints[11]);

  if (FMS_.FMSInit != 0) {
     throw std::runtime_error("unexpected FMS init flag");
  }

  if (mpiStatus.MPI_TAG != 2) {
     throw std::runtime_error("Only FMS mode is allowed in FMS mode");
  }

  if (FMS_.iState >= FMSNumStates) {
     throw std::runtime_error("Desired state higher than FMSNumStates");
  }

  // Work out which derivatives we need
  // DH for Surface Hopping in ABIN, TrajID==0
  if (FMS_.TrajID == 0) {
    printf("Receive derivative matrix logic from ABIN\n");
    recvCount = FMSNumStates * (FMSNumStates-1)/2 + FMSNumStates;
    MPI_Recv(bufints, recvCount, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
    checkRecvTag(mpiStatus);
    checkRecvCount(&mpiStatus, MPI_INT, recvCount);
    for (int i=0, ij=0; i<FMSNumStates; i++) {
      for (int j=i; j<FMSNumStates; j++, ij++) {
        FMS_.Derivs[ij] = bufints[ij];
        printf("%d ",FMS_.Derivs[ij] ? 1 : 0);
      }
      printf("\n");
    }
  }

  printf("Running trajectory %d on State %d\n", FMS_.TrajID, FMS_.iState);

  recvCount = 1;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus );
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  Data_.sim_time = bufdoubles[0];
  printf("Simulation time: %16.3f\n", Data_.sim_time);

  // Receive QM coordinates from ABIN in Bohrs
  printf("Receiving coordinates from ABIN\n");
  recvCount = FMS_.NAtoms * 3;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  for (int i = 0; i < recvCount; i++) {
    coordinates[i] = bufdoubles[i];
    printf("%g ", bufdoubles[i]);
  }
  printf("\n");

  printf("Receiving previous MOs\n");
  // Receiving previous diabatic MOs
  recvCount =  Data_.nbf * Data_.nbf;
  MPI_Recv(OldData_.MOs, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  if (FMS_.OldWfn == 1 && FMS_.iCall > 0) {
    check_array_equality(OldData_.MOs, Data_.MOs, Data_.nbf * Data_.nbf);
  }
  if (FMS_.OldWfn == 1 && FMS_.iCall == 0) {
    populate_array(Data_.MOs, Data_.nbf * Data_.nbf, 2.5, 2.0);
    check_array_equality(OldData_.MOs, Data_.MOs, Data_.nbf * Data_.nbf);
  }

  // Receiving previous CI vectors
  printf("Receiving previous CI vectors\n");
  recvCount =  FMSNumStates * Data_.nci;
  MPI_Recv(OldData_.CIvecs, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  if (FMS_.OldWfn == 1 && FMS_.iCall > 0) {
    check_array_equality(OldData_.CIvecs, Data_.CIvecs, FMSNumStates*Data_.nci);
  } else  if (FMS_.OldWfn == 1 && FMS_.iCall == 0) {
    populate_CIvecs(Data_.CIvecs, Data_.nci);
    check_array_equality(OldData_.CIvecs, Data_.CIvecs, FMSNumStates*Data_.nci);
  }

  printf("Receiving blob\n");
  recvCount = Data_.nblob;
  MPI_Recv(OldData_.Blob, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);
  if (FMS_.OldWfn == 1 && FMS_.iCall > 0) {
    check_array_equality(OldData_.Blob, Data_.Blob, Data_.nblob);
  } else  if (FMS_.OldWfn == 1 && FMS_.iCall == 0) {
    populate_array(Data_.Blob, Data_.nblob, 5.0, 3.0);
    check_array_equality(OldData_.Blob, Data_.Blob, Data_.nblob);
  }

  printf("Receiving real velocity vector for projected couplings\n");
  recvCount = FMS_.NAtoms * 3;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);

  printf("Receiving imaginary velocity vector for projected couplings\n");
  recvCount = FMS_.NAtoms * 3;
  MPI_Recv(bufdoubles, recvCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, abin_client, &mpiStatus);
  checkRecvTag(mpiStatus);
  checkRecvCount(&mpiStatus, MPI_DOUBLE, recvCount);

  FMS_.iCall++;

  return tag;
}


void TCServerMock::fms_send() {
  // TODO: This might not work for QM/MM
  int total_atoms = FMS_.NAtoms;
  printf("Sending Energies via MPI\n");
  MPI_Send(Data_.Energy, FMSNumStates, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);

  printf("Sending Transition Dipoles via MPI\n");
  MPI_Send(Data_.TDip, 3*(FMSNumStates-1), MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);

  printf("Sending Dipole Moment via MPI\n");
  MPI_Send(Data_.Dip, 3*FMSNumStates, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client); 

  printf("Sending Partial Charges via MPI\n");
  sendQMCharges();

  printf("Sending MOs via MPI.\n");
  MPI_Send(Data_.MOs, Data_.nbf*Data_.nbf, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);

  printf("Sending CI vectors via MPI.\n");
  MPI_Send(Data_.CIvecs, FMSNumStates*Data_.nci, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);

  printf("Sending wavefunction overlap via MPI.\n");
  MPI_Send(Data_.SMatrix, FMSNumStates*FMSNumStates, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);

  printf("Sending blob via MPI.\n");
  MPI_Send(Data_.Blob, Data_.nblob, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);

  printf("Sending derivative matrix via MPI\n");
  for (int i=0, ij=0; i<FMSNumStates; i++) {
    for (int j=i; j<FMSNumStates; j++, ij++) {
      int offset = ij * 3 * total_atoms;
      MPI_Send(&(Data_.DerivMat[offset]), 3*total_atoms, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
    }
  }

  // Currently not implemented in ABIN
  bool FMS_GAIMS = false;
  if (FMS_GAIMS) {
    printf("Sending SOMat via MPI\n");
    MPI_Send( Data_.SOMat, 18*FMSNumStates*FMSNumStates, MPI_DOUBLE, 0, MPI_TAG_OK, abin_client);
  }
}

fms_data::fms_data()
{
  // For now, let's just set these to a fixed values here.
  nbf = 10;
  nci = 9;
  nblob = 11;

  Energy = NULL;
  Dip = NULL;
  TDip = NULL;
  Chg = NULL;
  DerivMat = NULL;
  MOs = NULL;
  CIvecs = NULL;
  SMatrix = NULL;
  Blob = NULL;

  SOMat = NULL; // GAIMS
  sim_time = 0.0;
}

fms_data::~fms_data()
{
  TCTensor::freeDouble1D(Energy);
  TCTensor::freeDouble1D(Dip);
  TCTensor::freeDouble1D(TDip);
  TCTensor::freeDouble1D(Chg);
  TCTensor::freeDouble1D(DerivMat);
  TCTensor::freeDouble1D(MOs);
  TCTensor::freeDouble1D(CIvecs);
  TCTensor::freeDouble1D(SMatrix);
  TCTensor::freeDouble1D(Blob);

  TCTensor::freeDouble1D(SOMat); // GAIMS
}

fms_directive::fms_directive() {
  FMSInit = 0;
  NAtoms = 0;
  MMAtoms = 0;
  DoCoup = 0;
  TrajID = 0;
  Cent1 = 0;
  Cent2 = 0;
  StateID = 0;
  OldWfn = 0;
  iState = 0;
  jState = 0;
  FirstCall = 0;
  FMSRestart = 0;
  Derivs = NULL;
  // DH: This one is extra here to actually keep track
  // of how many times ABIN called us.
  iCall = 0;
}

fms_directive::~fms_directive() {
  free(Derivs);
}

void TCServerMock::allocate_fms_data() {
  if (FMSNumStates < 1) {
    throw std::runtime_error("FMSNumStates not set");
  }
  FMS_.Derivs = (bool*) malloc(sizeof(bool) * (FMSNumStates*(FMSNumStates+1)/2));
  memset(FMS_.Derivs, 0, sizeof(bool) * (FMSNumStates*(FMSNumStates+1)/2));

  Data_.Energy = TCTensor::double1D(FMSNumStates);
  Data_.Dip = TCTensor::double1D(FMSNumStates*3);
  Data_.TDip = TCTensor::double1D((FMSNumStates-1)*3);
  Data_.Chg = TCTensor::double1D(FMS_.NAtoms);
  Data_.DerivMat = TCTensor::double1D(FMSNumStates*(FMSNumStates+1) * 3 * FMS_.NAtoms / 2);

  Data_.MOs = TCTensor::double1D(Data_.nbf * Data_.nbf);
  Data_.CIvecs = TCTensor::double1D(FMSNumStates * Data_.nbf);
  Data_.SMatrix = TCTensor::double1D(FMSNumStates * FMSNumStates);
  Data_.Blob = TCTensor::double1D(Data_.nblob);

  // Currently not implemented in ABIN
  bool FMS_GAIMS = false;
  if (FMS_GAIMS)
    Data_.SOMat = TCTensor::double1D(18*FMSNumStates*FMSNumStates);
 
  OldData_.MOs = TCTensor::double1D(Data_.nbf * Data_.nbf);
  OldData_.CIvecs = TCTensor::double1D(FMSNumStates * Data_.nci);
  OldData_.SMatrix = TCTensor::double1D(FMSNumStates * FMSNumStates);
  OldData_.Blob = TCTensor::double1D(Data_.nblob);
}

void TCServerMock::populate_fms_data() {
  printf("Entering populate_fms_data\n");

  double energy = getWaterGradients();

  // Energies
  // For now just a constant shift between states
  double energyShift = 1.0; // in atomic units
  for (int i = 0; i < FMSNumStates; i++) {
    Data_.Energy[i] = energy + i*energyShift;
  }

  for (int i = 0; i < FMSNumStates * 3; i++) {
    Data_.Dip[i] = 0.0;
  }
  // For now populate only S0 dipole moment, reuse stuff from Amber interface
  double DTotal = 0.0;
  computeFakeQMDipoleMoment(Data_.Dip[0], Data_.Dip[1], Data_.Dip[2], DTotal);

  for (int i=0; i<FMSNumStates-1; i++) {
    // For now we only support singlets here.
    // if (FMS_Mults[i+1] == 1) {
    int FMS_Mults = 1;
    if (FMS_Mults == 1) {
        Data_.TDip[i*3+0] = i + 0.0;
        Data_.TDip[i*3+1] = i + 1.0;
        Data_.TDip[i*3+2] = i + 2.0;
    }
    // Assumes spin-pure states so singlet-triplet tdips are zero
    else {
        Data_.TDip[i*3+0] = 0.0;
        Data_.TDip[i*3+1] = 0.0;
        Data_.TDip[i*3+2] = 0.0;
    }
  }

  for (int i = 0, ij = 0; i < FMSNumStates; i++) {
    for (int j = 0; j < FMSNumStates; j++, ij++) {
      // Best case scenario for now
      // ABIN is ignoring this at the moment anyway
      if (i == j) {
        Data_.SMatrix[ij] = 1.0;
      } else {
        Data_.SMatrix[ij] = 0.0;
      }
    }
  }

  // TODO: Gradients/Couplings
  // For now just populate S0 gradient
  int ij = 0;
  int total_atoms = FMS_.NAtoms + FMS_.MMAtoms;
  int offset = ij * 3 * total_atoms;
  for(int i = 0; i < total_atoms; i++) {
    Data_.DerivMat[3*i + offset]   = gradients[3*i];
    Data_.DerivMat[3*i+1 + offset] = gradients[3*i+1];
    Data_.DerivMat[3*i+2 + offset] = gradients[3*i+2];
  }

  // TODO: It would be great if the data arrays were different in different time steps,
  // but it's tricky to get orking while also supporting restarting the simulation.
  populate_array(Data_.MOs, Data_.nbf * Data_.nbf, 2.5, 2.0);
  populate_array(Data_.Blob, Data_.nblob, 5.0, 3.0);
  populate_CIvecs(Data_.CIvecs, Data_.nci);
}

void TCServerMock::check_array_equality(double *A, double *B, int size) {
  for (int i = 0; i < size; i++) {
    if (A[i] != B[i]) {
      printf("%g != %g\n", A[i], B[i]);
      throw std::runtime_error("array validation failed");
    }
  }
}

void TCServerMock::print_array(double *A, int size) {
  for (int i = 0; i < size; i++) {
    printf("%g ", A[i]);
  }
  printf("\n");
}

void TCServerMock::populate_array(double *A, int size, double shift, double stride) {
  for (int i = 0; i < size; i++) {
    A[i] = shift + stride * i;
  }
}

void TCServerMock::populate_CIvecs(double *CIvecs, int nci) {
  for (int i = 0; i < FMSNumStates * nci; i++) {
    CIvecs[i] = 0.0;
  }
  CIvecs[0] = 1.0;
  CIvecs[nci + 1] = 1.0;
  CIvecs[2 * nci + 2] = 1.0;
}
