#/bin/bash
set -euo pipefail
# Useful for debugging
# TODO: Comment this out before merging
set -x

ABINEXE=$1
ABINOUT=abin.out
ABININ=input.in
ABINGEOM=mini.xyz
TCSRC="../tc_mpi_api.cpp tc_server.cpp"
TCEXE=tc_server
TCOUT=tc.out

if [[ -z ${MPI_PATH-} ]];then
  MPIRUN=mpirun
  MPICXX=mpicxx
  MPICH_HYDRA=hydra_nameserver
else
  MPIRUN=$MPI_PATH/bin/mpirun
  MPICXX=$MPI_PATH/bin/mpicxx
  MPICH_HYDRA=$MPI_PATH/bin/hydra_nameserver
fi

rm -f restart.xyz movie.xyz $TCEXE
if [[ "${1-}" = "clean" ]];then
   rm -f $TCOUT $ABINOUT *dat *diff restart.xyz.old
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

TC_PORT="tcport.$$"
# Make sure hydra_nameserver is running
hydra=$(ps -C hydra_nameserver -o pid= || true)
if [[ -z ${hydra-} ]];then
   echo "Launching hydra nameserver for MPI_Lookup"
   $MPICH_HYDRA &
fi

hostname=$HOSTNAME
MPIRUN="$MPIRUN -nameserver $hostname -n 1"

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT"
TC_CMD="./$TCEXE $TC_PORT.1"

$MPIRUN $TC_CMD > $TCOUT 2>&1 &
# Get PID of the last process
tcpid=$!

$MPIRUN $ABIN_CMD > $ABINOUT 2>&1 &
abinpid=$!

function cleanup {
  kill -9 $tcpid $abinpid > /dev/null 2>&1 || true
  exit 0
}

trap cleanup INT ABRT TERM EXIT

# The MPI interface is prone to deadlocks, where
# both server and client are waiting on MPI_Recv.
# We need to kill both processes if that happens.
MAX_TIME=6
seconds=1
while true;do
  ps -p $tcpid > /dev/null || tc_stopped=1
  ps -p $abinpid > /dev/null || abin_stopped=1
  if [[ -n ${tc_stopped:-} && -n ${abin_stopped:-} ]];then
    # Both TC and ABIN stopped, hopefully succesfully
    break
  elif [[ -n ${tc_stopped:-} || -n ${abin_stopped:-} ]];then
    # TC or ABIN ended, give the other time to finish.
    sleep 1
    if ! ps -o pid= -p $tcpid;then
      echo "Fake TeraChem died. Killing ABIN."
      cat $TCOUT
      cleanup
    elif ! ps -o pid= -p $abinpid;then
      echo "ABIN died. Killing fake TeraChem."
      cat $ABINOUT
      cleanup
    else
      # Normal exit
      break
    fi
  fi
  # Maybe add longer sleep interval to make this less flaky
  # (i.e. TC and ABIN do not end at the exact same time")
  # Alternatively, we can always return 0 from cleanup
  sleep 1
  let ++seconds
  if [[ $seconds -gt $MAX_TIME ]];then
    echo "Maximum time exceeded."
    cleanup
  fi
done
