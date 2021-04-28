#/bin/bash

# Test that ABIN exits gracefully
# when reading from empty TC port file.
touch port.txt.1

set -euo pipefail

source ../test_tc_server_utils.sh

IDX=7
ABININ=input.in$IDX
ABINOUT=${ABINOUT}$IDX

launch_hydra_nameserver $MPICH_HYDRA

MPIRUN="$MPIRUN -n 1"

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM"


$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 || true &
abinpid=$!

function cleanup {
  kill -9 $abinpid > /dev/null 2>&1 || true
  rm -f port.txt.1
  if [[ -f ERROR ]];then
    mv ERROR ABIN_ERROR$IDX
  fi
  exit 0
}

trap cleanup INT ABRT TERM EXIT

check_running_processes $abinpid
