#/bin/bash
set -euo pipefail

ABINEXE=$1
source ../test_tc_server_utils.sh

set_default_vars
set_mpich_vars
# If $1 = "clean"; exit early.
if ! clean_output_files $1; then
  exit 0
fi

N_TERA_SERVERS=$(egrep --only-matching 'nteraservers\s*=\s*[0-9]' $ABININ | egrep -o [0-9])

# Exit early for OpenMPI build.
check_for_openmpi

# Compile fake TC server
$MPICXX $TCSRC -Wall -o $TCEXE

# NOTE: We very intentionally do NOT launch
# hydra_nameserver in this test since it cannot handle
# multiple TC servers due to a bug in MPI_Unpublish_name
# https://github.com/pmodels/mpich/issues/5058
#
# Therefore, we pass the port_name to ABIN via files, see below.

MPIRUN="$MPIRUN -n 1"

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM"

function cleanup {
  kill -9 ${job_pids[@]} > /dev/null 2>&1 || true
  exit 0
}
trap cleanup INT ABRT TERM EXIT

let NUM_JOBS=N_TERA_SERVERS+1
declare -A job_pids
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
   $MPIRUN ./$TCEXE > $TCOUT.$itera 2>&1 &
   job_pids[$itera]=$!
}
sleep 2
sync

# Grep port names from TC output, pass to ABIN via a file.
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
  grep 'port name' $TCOUT.$itera | awk -F"port name: " '{print $2;exit}' > $TC_PORT_FILE.$itera
}

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 &
job_pids[$NUM_JOBS]=$!

check_running_processes ${job_pids[@]}
