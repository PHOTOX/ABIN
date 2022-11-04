#/bin/bash
set -euo pipefail

export ABINEXE=$1

source ../test_tc_server_utils.sh

set_default_vars
set_mpich_vars

if ! clean_output_files $1; then
  exit 0
fi

check_for_openmpi

# Compile default TC server
$MPICXX $TCSRC -Wall -o $TCEXE

launch_hydra_nameserver $MPICH_HYDRA

MPIRUN="$MPIRUN -nameserver $HOSTNAME -n 1"

TC_PORT="tcport.$$"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT -v veloc.in"
TC_CMD="./$TCEXE $TC_PORT.1"

function cleanup {
  kill -9 $tcpid $abinpid > /dev/null 2>&1 || true
  exit 0
}
trap cleanup INT ABRT TERM EXIT

$MPIRUN $TC_CMD > $TCOUT 2>&1 &
tcpid=$!

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 &
abinpid=$!

check_running_processes $abinpid $tcpid

