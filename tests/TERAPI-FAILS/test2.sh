#/bin/bash

# Test that ABIN sends TC error tag early
# when it fails during parsing its own input,
# so that the TC server exits gracefully.

set -euo pipefail

source ../test_tc_server_utils.sh

ABININ=input.in2
ABINOUT=${ABINOUT}2
TCOUT=${TCOUT}2

launch_hydra_nameserver $MPICH_HYDRA
hostname=$HOSTNAME
MPIRUN="$MPIRUN -nameserver $hostname -n 1"

TC_PORT="test2.$$"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT"
TC_CMD="./$TCEXE $TC_PORT.1"

$MPIRUN $TC_CMD > $TCOUT 2>&1 || true &
tcpid=$!

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 || true &
abinpid=$!

function cleanup {
  kill -9 $tcpid $abinpid > /dev/null 2>&1 || true
  grep 'what()' $TCOUT > TC_ERROR2
  exit 0
}

trap cleanup INT ABRT TERM EXIT

check_running_processes $abinpid $tcpid
