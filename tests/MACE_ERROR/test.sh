#!/bin/bash
set -uo pipefail

ABINEXE=$1

ABINOUT=abin.out
ABININ=input.in
ABINGEOM=mini.xyz

MACE_SERVER=../../interfaces/MACE/mace_server.py
MACE_OUT=mace_server.out

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM"
MACE_CMD="python3 $MACE_SERVER --device cpu --model-path __MOCK_ERROR__"
PORT_FILE=mace_port.txt

# If $1 = "clean"; exit early.
if [[ "${1-}" = "clean" ]]; then
  rm -f $MACE_OUT $ABINOUT ./*.dat ./*.diff
  rm -f restart.xyz velocities.xyz forces.xyz movie.xyz restart.xyz.old
  rm -f $PORT_FILE ERROR ompi_uri.txt
  exit 0
fi

# Skip the test if the python environment does not contain necessary libraries
if ! python3 -c "import mpi4py, numpy, ase" 2>/dev/null; then
    echo "MACE python environment not activated"
    exit 3
fi

# Determine MPI paths
if [[ -z ${MPI_PATH-} ]]; then
  MPIRUN=mpirun
else
  MPIRUN=$MPI_PATH/bin/mpirun
  export PATH=$MPI_PATH/bin:$PATH
  export LD_LIBRARY_PATH=$MPI_PATH/lib:$LD_LIBRARY_PATH
fi

# Detect OpenMPI vs MPICH
IS_OPENMPI=false
if $MPIRUN --version 2>&1 | grep -q "Open MPI"; then
  IS_OPENMPI=true
fi

ompi_server_pid=""
MPIRUN_EXTRA_ARGS=""

if [[ "$IS_OPENMPI" = "true" ]]; then
  # OpenMPI requires ompi-server for MPI_Comm_connect/accept
  OMPI_SERVER=${MPI_PATH-}/bin/ompi-server
  if [[ -z ${MPI_PATH-} ]]; then
    OMPI_SERVER=ompi-server
  fi

  if ! which $OMPI_SERVER &> /dev/null; then
    echo "ERROR: Skipping MACE test: ompi-server not found (required for OpenMPI)" >> ERROR
    exit 1
  fi

  OMPI_URI_FILE="$PWD/ompi_uri.txt"
  $OMPI_SERVER --no-daemonize -r "$OMPI_URI_FILE" &
  ompi_server_pid=$!
  sleep 1

  if [[ ! -f "$OMPI_URI_FILE" ]]; then
    echo "ERROR: ompi-server did not create URI file" >> ERROR
    kill $ompi_server_pid 2> /dev/null
    exit 1
  fi

  MPIRUN_EXTRA_ARGS="--ompi-server file:$OMPI_URI_FILE"
fi

MPIRUN_CMD="$MPIRUN -n 1 $MPIRUN_EXTRA_ARGS"

# Cleanup function to stop the background processes
function cleanup {
  if [[ -n ${mace_pid-} ]] && kill -0 $mace_pid >& /dev/null; then
    echo "ERROR: MACE server $mace_pid is still running!" >> ERROR
    kill ${mace_pid-} &> /dev/null || true
  fi

  if [[ -n ${abinpid-} ]] && kill -0 $abinpid >& /dev/null; then
    echo "ERROR: ABIN process $abinpid is still running!" >> ERROR
    kill ${abinpid-} &> /dev/null || true
  fi

  if [[ -n "${ompi_server_pid-}" ]]; then
    kill $ompi_server_pid &> /dev/null || true
  fi
}

function wait_for_portfile {
  # Wait 10s for the MACE server to write the port file
  MAX_WAIT=20
  i=0
  while [[ ! -f $PORT_FILE ]]; do
    if [[ $i -gt $MAX_WAIT ]] || ! kill -0 $mace_pid >& /dev/null; then
      echo "ERROR: MACE server did not write port file 'mace_port.txt'" | tee ERROR
      set -x
      cat $MACE_OUTPUT
      exit 1
    fi
    sleep 0.5
    let i++
  done
}

# Automatically call the cleanup function when the script
# exits or is interrupted by a signal
trap cleanup INT ABRT TERM EXIT

# Launch mock MACE server
$MPIRUN_CMD $MACE_CMD &> $MACE_OUT &
mace_pid=$!

wait_for_portfile

# Launch ABIN
$MPIRUN_CMD $ABIN_CMD &> $ABINOUT &
abinpid=$!

# Give both processes 10 seconds to finish
MAX_ITER=20
iter=0
# Note about 'kill -0' https://unix.stackexchange.com/questions/169898/what-does-kill-0-do
while ( (kill -0 $abinpid >& /dev/null) || (kill -0 $mace_pid >& /dev/null) ); do
  if [[ $iter -gt $MAX_ITER ]]; then
    echo "Test did not finish in time" >> ERROR
    break
  fi
  sleep 0.5
  let iter++
done
