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


# TODO: Simplify this
tc_stopped=
abin_stopped=
function check_tc {
  tc_stopped=
  ps -p $tcpid > /dev/null || tc_stopped=1
}

function check_abin {
  abin_stopped=
  ps -p $abinpid > /dev/null || abin_stopped=1
}


# TODO: Move this to ../tc_
# The MPI interface is prone to deadlocks, where
# both server and client are waiting on MPI_Recv.
# We need to kill both processes if that happens.
MAX_ITER=100
iter=1
while true;do
  check_abin
  check_tc
  if [[ -n ${tc_stopped:-} && -n ${abin_stopped:-} ]];then
    # Both TC and ABIN stopped.
    break
  elif [[ -n ${tc_stopped:-} || -n ${abin_stopped:-} ]];then
    # TC or ABIN ended, give the other time to finish.
    sleep 1
    check_abin
    check_tc
    if [[ -n ${tc_stopped:-} && -n ${abin_stopped:-} ]];then
      # Both TC and ABIN stopped.
      break
    elif [[ -n ${tc_stopped:-} ]];then
      echo "Fake TeraChem died. Killing ABIN."
      echo "Printing TC and ABIN outputs"
      cat $TCOUT $ABINOUT
      cleanup
    else
      echo "ABIN died. Killing fake TeraChem."
      echo "Printing TC and ABIN outputs"
      cat $TCOUT $ABINOUT
      cleanup
    fi
  fi

  sleep 0.2
  let ++iter
  if [[ $iter -gt $MAX_ITER ]];then
    echo "Maximum time exceeded."
    cleanup
  fi
done
