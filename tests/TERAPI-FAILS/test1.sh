#/bin/bash

# Test a scenario in which SCF does not converge
# and TeraChem sends the MPI_SCF_DIE tag to ABIN.
# In that case, ABIN should stop with an error.

set -euo pipefail

source ../test_tc_server_utils.sh

IDX=1
ABININ=input.in$IDX
ABINOUT=${ABINOUT}$IDX
TCOUT=${TCOUT}$IDX
TCSRC="../tc_mpi_api.cpp ../../water_potentials/qtip4pf.cpp tc_server$IDX.cpp"
TCEXE=tc_server$IDX

# Compile fake TC server
$MPICXX $TCSRC -Wall -o $TCEXE

launch_hydra_nameserver $MPICH_HYDRA

hostname=$HOSTNAME
MPIRUN="$MPIRUN -nameserver $hostname -n 1"

TC_PORT="test1.$$"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT"
TC_CMD="./$TCEXE $TC_PORT.1"

$MPIRUN $TC_CMD > $TCOUT 2>&1 || true &
tcpid=$!

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 || true &
abinpid=$!

function cleanup {
  kill -9 $tcpid $abinpid > /dev/null 2>&1 || true
  grep 'what()' $TCOUT > TC_ERROR$IDX
  if [[ -f ERROR ]];then
    mv ERROR ABIN_ERROR$IDX
  fi
  exit 0
}

trap cleanup INT ABRT TERM EXIT

check_running_processes $abinpid $tcpid
