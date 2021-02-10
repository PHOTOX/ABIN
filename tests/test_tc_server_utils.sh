#!/bin/bash

# Various utility function that are used by 
# tests for TeraChem MPI interface (e.g. in TERAPI/)

# This file is meant to be sourced, NOT executed!

set -euo pipefail

export TC_PORT_FILE=port.txt

hydrapid=
launch_hydra_nameserver() {
  # Make sure hydra_nameserver is running
  CMD=$1
  hydra=$(ps -C hydra_nameserver -o pid= || true)
  if [[ -z ${hydra-} ]];then
    #echo "Launching hydra nameserver for MPI_Lookup"
    $CMD &
    hydrapid=$!
  fi
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
  rm -f $TCEXE $TCOUT* $ABINOUT $TC_PORT_FILE.* ERROR
  return $return_code
}


# Sillently kill all processes whose PIDs
# are passed as parameters. Note that some of
# them could have already ended sucessfully.
kill_processes() {
  kill -9 $* > /dev/null 2>&1 || true
}

# Not that it is hard in general to know
# whether we ended successfully or not,
# so we always return 0. Validation is then
# always done on the output files.
cleanup() {
  kill_processeses $*
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
    running=$(ps -eo pid|grep -E "$regex"|wc -l)
    if [[ $running -eq 0 ]];then
      #echo "Both ABIN and TeraChem servers stopped"
      break
    elif [[ $running -lt $num_jobs ]];then
      # Give the others time to finish 
      sleep 1
      running=$(ps -eo pid|grep -E "$regex"|wc -l)
      if [[ $running -ne 0 ]];then
        echo "One of the TC servers or ABIN died. Killing the rest."
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
