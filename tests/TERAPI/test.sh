#/bin/bash
set -euo pipefail
# Useful for debugging
#set -x

ABINEXE=$1
source ../test_tc_server_utils.sh

set_default_vars
set_mpich_vars
# If $1 = "clean"; exit early.
if ! clean_output_files $1; then
  exit 0
fi

# Exit early for OpenMPI build.
check_for_openmpi

# Compiled the fake TC server
$MPICXX $TCSRC -Wall -o $TCEXE

launch_hydra_nameserver $MPICH_HYDRA

hostname=$HOSTNAME
MPIRUN="$MPIRUN -nameserver $hostname -n 1"

TC_PORT="tcport.$$"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT"
TC_CMD="./$TCEXE $TC_PORT.1"

$MPIRUN $TC_CMD > $TCOUT 2>&1 &
tcpid=$!

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 &
abinpid=$!

function cleanup {
  kill -9 $tcpid $abinpid > /dev/null 2>&1 || true
  #kill -9 $tcpid $abinpid $hydrapid > /dev/null 2>&1 || true
  exit 0
}

trap cleanup INT ABRT TERM EXIT

check_running_processes $abinpid $tcpid
