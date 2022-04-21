#/bin/bash
set -euo pipefail

function local_cleanup {
  kill -9 ${job_pids[@]}  > /dev/null 2>&1 || true
  exit 0
}

ABINEXE=$1
source ../test_tc_server_utils.sh

set_default_vars
set_mpich_vars
rm -f remd.out abin.out restart.xyz.??.old restart.xyz.?? restart.xyz.??.? geom.dat.?? movie.xyz.?? temper.dat.?? energies.dat.?? 
if ! clean_output_files $1; then
  exit 0
fi

N_TERA_SERVERS=$(egrep --only-matching 'nreplica\s*=\s*[0-9]' $ABININ | egrep -o [0-9])

# Exit early for OpenMPI build.
check_for_openmpi

# Compiled the fake TC server
$MPICXX $TCSRC -Wall -o $TCEXE

# NOTE: We very intentionally do NOT launch
# hydra_nameserver in this test since it cannot handle
# multiple TC servers due to a bug in MPI_Unpublish_name
# https://github.com/pmodels/mpich/issues/5058
TC_SERVER_NAME="tcserver.$$"

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM" # -M $TC_SERVER_NAME"

trap local_cleanup INT ABRT TERM EXIT

let NUM_JOBS=N_TERA_SERVERS+1
declare -A job_pids
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
   $MPIRUN -n 1 ./$TCEXE $TC_SERVER_NAME.$itera > $TCOUT.$itera 2>&1 &
   job_pids[$itera]=$!
}
sleep 2
# Grep port names from TC output, pass to ABIN via a file.
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
  grep 'port name' $TCOUT.$itera | awk -F"port name: " '{print $2;exit}' > $TC_PORT_FILE.$itera
}

$MPIRUN -np 3 $ABINEXE -i input.in -x mini.xyz -v vel0.in > $ABINOUT 2>&1
job_pids[$NUM_JOBS]=$!

check_running_processes ${job_pids[@]}
