#/bin/bash

# Test that ABIN stops if nteraservers > nbeads

set -euo pipefail
source ../test_tc_server_utils.sh

IDX=3
ABININ=input.in$IDX
ABINOUT=${ABINOUT}$IDX
TCOUT=${TCOUT}$IDX
N_TERA_SERVERS=$(egrep --only-matching 'nteraservers\s*=\s*[0-9]' $ABININ | egrep -o [0-9])

MPIRUN="$MPIRUN -n 1"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM"

declare -A job_pids
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
   $MPIRUN ./$TCEXE > $TCOUT.$itera 2>&1 &
   job_pids[$itera]=$!
}
sleep 1
# Grep port names from TC outputs, pass to ABIN via a file.
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
  grep 'port name' $TCOUT.$itera | awk -F"port name: " '{print $2;exit}' > $TC_PORT_FILE.$itera
}

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 &
job_pids[$(expr $N_TERA_SERVERS + 1)]=$!

function cleanup {
  kill -9 ${job_pids[@]} > /dev/null 2>&1 || true
  grep 'what()' $TCOUT.* > TC_ERROR$IDX
  if [[ -f ERROR ]];then
    mv ERROR ABIN_ERROR$IDX
  fi
  exit 0
}

trap cleanup INT ABRT TERM EXIT

check_running_processes ${job_pids[@]}
