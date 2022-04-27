#!/bin/bash

# Various utility function that are shared by test.sh scripts
# in tests for TeraChem MPI interface (TERAPI*/)

# This file is meant to be sourced, NOT executed!

set -euo pipefail

export TC_PORT_FILE=port.txt

hydrapid=
launch_hydra_nameserver() {
  # Make sure hydra_nameserver is running
  # NOTE: Currently we will always restart the server
  # to workaround the existing bug in it.
  # https://github.com/pmodels/mpich/issues/5058
  CMD=$1
  hydra=$(ps -C hydra_nameserve -o pid= || true)
  if [[ -n ${hydra} ]];then
    kill_processes $hydra
  fi
  #if [[ -z ${hydra-} ]];then
    #echo "Launching hydra nameserver for MPI_Lookup"
  #fi
  $CMD &
  hydrapid=$!
  # Sometime tests fail connecting to the nameserver,
  # let's try to give it some time.
  sleep 0.2
}

check_for_openmpi() {
  if [[ -f "${MPI_PATH-}/bin/orterun" ]];then
    # TeraChem is compiled with MPICH so there's no
    # point in trying to make this work with OpenMPI.
    # We'll skip this test by faking it was successfull.
    # Here are some pointers if we ever want to make it work:
    # https://techdiagnosys.blogspot.com/2016/12/openmpi-working-nameserver-publish.html
    # https://www.open-mpi.org/doc/v4.1/man1/ompi-server.1.php
    # https://www.open-mpi.org/doc/v4.1/man1/mpirun.1.php#sect6 (search for ompi-server)
    echo "Skipping this test with OpenMPI build"
    for f in `ls *ref`;do
      cp $f `basename $f .ref`
    done
    exit 1
  fi
}

check_for_intelmpi() {
  if which mpiifort > /dev/null;then
    echo "Skipping this test for IntelMPI build"
    for f in `ls *ref`;do
      cp $f `basename $f .ref`
    done
    exit 1
  fi
}

set_mpich_vars() {
  if [[ -z ${MPI_PATH-} ]];then
    export MPIRUN=mpirun
    export MPICXX=mpicxx
    export MPICH_HYDRA=hydra_nameserver
  else
    export MPIRUN=$MPI_PATH/bin/mpirun
    export MPICXX=$MPI_PATH/bin/mpicxx
    export MPICH_HYDRA=$MPI_PATH/bin/hydra_nameserver
  fi
}

set_default_vars() {
  export ABINOUT=abin.out
  export ABININ=input.in
  export ABINGEOM=mini.xyz
  export TCSRC="../tc_mpi_api.cpp ../../water_potentials/qtip4pf.cpp tc_server.cpp"
  export TCEXE=tc_server
  export TCOUT=tc.out
}

clean_output_files() {
  local return_code=0
  if [[ -z "${1-}" ]];then
    echo "Incorrect invocation of the clean_output_files function!"
    return 1
  fi
  if [[ "${1-}" = "clean" ]];then
    return_code=1
  fi
  # WARNING: This shift is important, since by default
  # the first parameter is the path to ABIN binary,
  # we do not want to delete that!
  shift
  # Remove predefined files + any additional arguments that were passed.
  rm -f $* 
  rm -f *dat *diff
  rm -f restart.xyz velocities.xyz forces.xyz movie.xyz restart.xyz.old
  rm -f $TCEXE $TCOUT* $ABINOUT $TC_PORT_FILE.* ERROR scrdir000?
  return $return_code
}


# Sillently kill all processes whose PIDs are passed as parameters.
# Typically, all these processes should have already ended.
kill_processes() {
  kill -9 $* > /dev/null 2>&1 || true
}

# NOTE that this function will typically get overwritten
# in the TERAPI*/test.sh scripts.
cleanup() {
  kill_processes $*
  # It is hard in general to know whether we ended successfully or not,
  # so we always return 0. Validation is then
  # always done on the output files.
  exit 0
}

# Helper function for building a regex expression
join_by() {
  local IFS="$1"
  shift
  echo "$*"
}

# This function accepts PIDs of ABIN and all TC servers
# end periodically checks whether they are still running.
# If only some of them stopped, it kills the rest.
check_running_processes() {
  # The MPI interface is prone to deadlocks, where
  # both server and client are waiting on MPI_Recv.
  # We need to kill both processes if that happens.
  pids="$*"
  num_jobs=$#
  regex=$(join_by \| $pids)
  MAX_ITER=100
  iter=1
  while true;do
    running=$(ps -eo pid | egrep "$regex" | wc -l)
    if [[ $running -eq 0 ]];then
      break
    elif [[ $running -lt $num_jobs ]];then
      # Give the others time to finish 
      sleep 1.2
      running=$(ps -eo pid | egrep "$regex" | wc -l)
      if [[ $running -ne 0 ]];then
        echo "One of the TC servers or ABIN died. Killing the rest."
	#echo "Printing ABIN and TC server outputs for debugging."
	#echo "##################################################"
        #cat ${ABINOUT}* ${TCOUT}*
      fi
      break
    fi
 
    sleep 0.2
    let ++iter
    if [[ $iter -gt $MAX_ITER ]];then
      echo "Maximum time exceeded."
      break
    fi
  done
  # ALWAYS call cleanup
  cleanup $pids
}
