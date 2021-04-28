#/bin/bash

# Test how ABIN handles MPI error.
# Concretely, how it handles if TC sends less data
# than expected (see tc_server8.cpp).

set -euo pipefail

source ../test_tc_server_utils.sh

IDX=8
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

TC_PORT="test$IDX.$$"
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
