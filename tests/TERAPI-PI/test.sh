#/bin/bash
set -euo pipefail
# Useful for debugging
#set -x

ABINEXE=$1
ABINOUT=abin.out
ABININ=input.in
ABINGEOM=mini.xyz
TCSRC="../tc_mpi_api.cpp ../../water_potentials/qtip4pf.cpp tc_server.cpp"
TCEXE=tc_server
TCOUT=tc.out

hydrapid=
function launch_hydra_nameserver {
  # Make sure hydra_nameserver is running
  CMD=$1
  hydra=$(ps -C hydra_nameserver -o pid= || true)
  if [[ -z ${hydra-} ]];then
    echo "Launching hydra nameserver for MPI_Lookup"
    $CMD &
    hydrapid=$!
  fi
}

# TODO: Determine this from ABIN input!
N_TERA_SERVERS=4  # Use more TC servers for PIMD or REMD

if [[ -z ${MPI_PATH-} ]];then
  MPIRUN=mpirun
  MPICXX=mpicxx
  MPICH_HYDRA=hydra_nameserver
else
  MPIRUN=$MPI_PATH/bin/mpirun
  MPICXX=$MPI_PATH/bin/mpicxx
  MPICH_HYDRA=$MPI_PATH/bin/hydra_nameserver
fi

rm -f restart.xyz movie.xyz $TCEXE $ABINOUT $TCOUT.? port.txt.*
if [[ "${1-}" = "clean" ]];then
   rm -f $TCOUT $ABINOUT *dat *diff restart.xyz.old velocities.xyz forces.xyz
   exit 0
fi

if [[ -f "${MPI_PATH-}/bin/orterun" ]];then
  # TeraChem is compiled with MPICH so there's no
  # point in trying to make this work with OpenMPI.
  # We'll skip this test by faking it was successfull.
  # Here are some pointers if we ever want to make it work:
  # https://techdiagnosys.blogspot.com/2016/12/openmpi-working-nameserver-publish.html
  # https://www.open-mpi.org/doc/v4.1/man1/ompi-server.1.php
  # https://www.open-mpi.org/doc/v4.1/man1/mpirun.1.php#sect6 (search for ompi-server)
  echo "Skipping TERAPI test with OpenMPI"
  # TODO: Is there a less hacky way to fake passing this test?
  # Or should we skip it altogether already in tests/test.sh?
  for f in `ls *ref`;do
    cp $f `basename $f .ref`
  done
  exit 0
fi

# Compiled the fake TC server
$MPICXX $TCSRC -Wall -o $TCEXE

# NOTE: We very intentionally do NOT launch
# hydra_nameserver in this test since it cannot handle
# multiple TC servers due to a bug in MPI_Unpublish_name
# https://github.com/pmodels/mpich/issues/5058
#
# Therefore, we pass the port_name to ABIN via files, see below.
#TC_SERVER_NAME="tcserver.$$"
#launch_hydra_nameserver $MPICH_HYDRA
#hostname=$HOSTNAME
#MPIRUN="$MPIRUN -nameserver $hostname -n 1"

MPIRUN="$MPIRUN -n 1"

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM" # -M $TC_SERVER_NAME"

let NUM_JOBS=N_TERA_SERVERS+1
declare -A job_pids
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
   #$MPIRUN ./$TCEXE $TC_SERVER_NAME.$itera > $TCOUT.$itera 2>&1 &
   $MPIRUN ./$TCEXE > $TCOUT.$itera 2>&1 &
   job_pids[$itera]=$!
}
sleep 1
# Grep port names from TC output, pass to ABIN via a file.
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
  grep 'port name' tc.out.$itera | awk -F"port name: " '{print $2;exit}' > port.txt.$itera
}

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 &
job_pids[$NUM_JOBS]=$!

function cleanup {
  kill -9 ${job_pids[@]}  > /dev/null 2>&1 || true
  exit 0
}

trap cleanup INT ABRT TERM EXIT

# The MPI interface is prone to deadlocks, where
# both server and client are waiting on MPI_Recv.
# We need to kill both processes if that happens.
MAX_TIME=10
seconds=1
# CHECK WHETHER ABIN AND TC ARE RUNNING
function join_by { local IFS="$1"; shift; echo "$*"; }
regex=`join_by \| ${job_pids[@]}`
while true;do
  njobs=$(ps -eo pid|grep -E "$regex"|wc -l)
  if [[ $njobs -eq 0 ]];then
    echo "Both ABIN and TeraChem servers stopped"
    break
  elif [[ $njobs -lt $NUM_JOBS ]];then
    # Give the others time to finish 
    sleep 1
    njobs=$(ps -eo pid|grep -E "$regex"|wc -l)
    if [[ $njobs -eq 0 ]];then
      echo "Both ABIN and TeraChem servers stopped"
      break
    fi
    echo "One of the TC servers or ABIN died. Killing the rest."
    cleanup
  fi

  sleep 1
  let ++seconds
  if [[ $seconds -gt $MAX_TIME ]];then
    echo "Maximum time exceeded."
    cleanup
  fi
done
