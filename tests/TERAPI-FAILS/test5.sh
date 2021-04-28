#/bin/bash

# Test that ABIN timeouts gracefully
# if it cannot find TC port via hydra_nameserver.

set -euo pipefail

source ../test_tc_server_utils.sh

IDX=5
ABININ=input.in$IDX
ABINOUT=${ABINOUT}$IDX

launch_hydra_nameserver $MPICH_HYDRA

hostname=$HOSTNAME
#MPIRUN="$MPIRUN -nameserver $hostname -n 1"

TC_PORT="test$IDX.$$"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT"

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 || true &
abinpid=$!

function cleanup {
  kill -9 $abinpid > /dev/null 2>&1 || true
  if [[ -f ERROR ]];then
    mv ERROR ABIN_ERROR$IDX
  fi
  exit 0
}

trap cleanup INT ABRT TERM EXIT

check_running_processes $abinpid
