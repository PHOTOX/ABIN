#/bin/bash

# Test that ABIN sends TC error tag early
# when it fails during parsing its own input,
# so that the TC server exits gracefully.

set -euo pipefail

source ../test_tc_server_utils.sh

function cleanup {
  kill -9 $tcpid $abinpid > /dev/null 2>&1 || true
  grep 'what()' $TCOUT > TC_ERROR$IDX
  if [[ -f ERROR ]];then
    mv ERROR ABIN_ERROR$IDX
  fi
  exit 0
}

IDX=2
ABININ=input.in$IDX
ABINOUT=${ABINOUT}$IDX
TCOUT=${TCOUT}$IDX

launch_hydra_nameserver $MPICH_HYDRA

MPIRUN="$MPIRUN -nameserver $HOSTNAME -n 1"

TC_PORT="test$IDX.$$"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT"
TC_CMD="./$TCEXE $TC_PORT.1"

trap cleanup INT ABRT TERM EXIT

$MPIRUN $TC_CMD > $TCOUT 2>&1 || true &
tcpid=$!

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 || true &
abinpid=$!

check_running_processes $abinpid $tcpid
