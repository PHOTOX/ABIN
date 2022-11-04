#include <cstdlib>
#include <math.h>
#include <cstring>
#include <stdexcept>
#include "mpi.h"

#include "../water_potentials/qtip4pf.h"

#define MAX_DATA 200

#define MPI_TAG_EXIT 0
#define MPI_TAG_GRADIENT 2
#define MPI_TAG_ERROR 13

// This is sent to ABIN for normal MPI_Send
// it is confusing that this is not symmetric, i.e.
// 0 from TC means OK, while 0 from ABIN means exit.
#define MPI_TAG_OK 0

#define MPI_TAG_SCF_DIED 1

// What does FMS want us to do?
// Store all the stuff sent by FMS
struct fms_directive {
  int FMSInit;
  int NAtoms;
  int MMAtoms;
  int DoCoup;
  int TrajID;
  int Cent1;
  int Cent2;
  int StateID;
  int OldWfn;
  int iState;
  int jState;
  int FirstCall;
  int FMSRestart;
  bool* Derivs;
  int iCall;

  fms_directive();
  ~fms_directive();
};

// All data that FMS/ABIN wants
struct fms_data {
  // Sizing
  // DH: Not sure why natoms is here, it's also in fms_directive
  // int natoms;
  int nbf;
  int nci;
  int nblob;

  // State information
  double *Energy;
  double *Dip;
  double *TDip;
  double *Chg;
  double *DerivMat;

  // Electronic structure info
  double *MOs;
  double *CIvecs;
  double *SMatrix;
  double *Blob;
  double sim_time;

  // GAIMS
  double* SOMat;

  fms_data();
  ~fms_data();
};

class TCServerMock {
  public:
    TCServerMock(char *);
    ~TCServerMock(void);

    void initializeCommunication();

    int receiveNumAtoms();
    char* receiveAtomTypes();
    void receiveAtomTypesAndScrdir();
    void receiveCoordinates();

    int receiveBeginLoop();
    // This combines all expected receive
    // methods in their order. But the tc_server
    // implementation can also call them individually
    // for better granurality.
    int receive();

    // Using qTIP4PF, same as in ABIN and used in all other tests.
    // Returns potential energy.
    double getWaterGradients();

    void send();
    void sendSCFEnergy(double, int);
    void sendQMCharges();
    void sendQMDipoleMoment();
    void sendQMGradients();

    // FMS Interface
    int fmsinit_receive();
    void fmsinit_send();
    int fms_receive();
    void fms_send();

    void allocate_fms_data();
    void populate_fms_data();
    void populate_fms_data_from_file(int);

    // Unfortunately, the number of FMS/SH states are given in TC input,
    // NOT via MPI interface. So we have no way to pass this from ABIN directly.
    // As a workaround, this value needs to be set manually in tc_server.cpp
    int FMSNumStates;

    MPI_Comm* getABINCommunicator();

  private:
    // This one will be published via MPI_Publish
    // ABIN will be looking for this one via hydra_nameserver
    char *tcServerName;
    // This is the name of the actual port from MPI_Open_port()
    char mpiPortName[MPI_MAX_PORT_NAME];
    // Buffers for MPI_Recv and MPI_Send calls
    char   bufchars[MAX_DATA];
    double bufdoubles[MAX_DATA];
    int    bufints[MAX_DATA];

    double *gradients;
    double *coordinates;

    // FMS / Surface Hopping stuff
    fms_data Data_;
    fms_data OldData_;
    fms_directive FMS_;

    MPI_Comm abin_client;
    MPI_Status mpiStatus;

    int totNumAtoms;
    char *atomTypes;

    void computeFakeQMDipoleMoment(double&, double&, double&, double&);
    void computeFakeQMCharges(double*);

    void publishServerName(char*, char*);
    void printMPIError(int);
    void parseScrDir(char*, char*);
    void validateScrDir(char*);
    void checkRecvCount(MPI_Status*, MPI_Datatype, int);
    void checkRecvTag(MPI_Status&);

    void populate_array(double*, int, double, double);
    void populate_CIvecs(double*, int);
    void check_array_equality(double*, double*, int);
    void print_array(double*, int);
};
