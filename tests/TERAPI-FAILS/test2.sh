#/bin/bash
set -euo pipefail
# Useful for debugging
#set -x

ABININ=input.in.wrong
TCSRC="../tc_mpi_api.cpp ../../water_potentials/qtip4pf.cpp tc_server1.cpp"

# Compiled fake TC server
$MPICXX $TCSRC -Wall -o $TCEXE

hostname=$HOSTNAME
MPIRUN="$MPIRUN -nameserver $hostname -n 1"

TC_PORT="test1.$$"
ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM -M $TC_PORT"
TC_CMD="./$TCEXE $TC_PORT.1"

$MPIRUN $TC_CMD > $TCOUT 2>&1 &
tcpid=$!

$MPIRUN $ABIN_CMD >> $ABINOUT 2>&1 &
abinpid=$!

function cleanup {
  kill -9 $tcpid $abinpid > /dev/null 2>&1 || true
  grep_tc_error $TCOUT
  exit 0
}

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

trap cleanup INT ABRT TERM EXIT

# The MPI interface is prone to deadlocks, where
# both server and client are waiting on MPI_Recv.
# We need to kill both processes if that happens.
MAX_TIME=6
seconds=1
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
    else
      cleanup
    fi
  fi

  sleep 1
  let ++seconds
  if [[ $seconds -gt $MAX_TIME ]];then
    echo "Maximum time exceeded."
    cleanup
  fi
done
