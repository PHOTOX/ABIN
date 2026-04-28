#!/bin/bash
set -euo pipefail

ABINEXE=$1

export ABINOUT=abin.out
export ABININ=input.in
export ABINGEOM=mini.xyz

MACE_SERVER=mock_mace_server.py
MACE_OUT=mace_server.out

# If $1 = "clean"; exit early.
if [[ "${1-}" = "clean" ]]; then
  rm -f $MACE_OUT $ABINOUT *.dat *.diff
  rm -f restart.xyz velocities.xyz forces.xyz movie.xyz restart.xyz.old
  rm -f mace_port.txt.* ERROR ompi_uri.txt
  exit 0
fi

# Determine MPI paths
if [[ -z ${MPI_PATH-} ]]; then
  export MPIRUN=mpirun
else
  export MPIRUN=$MPI_PATH/bin/mpirun
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

  if ! which $OMPI_SERVER > /dev/null 2>&1; then
    echo "Skipping MACE test: ompi-server not found (required for OpenMPI)"
    exit 0
  fi

  OMPI_URI_FILE="$PWD/ompi_uri.txt"
  $OMPI_SERVER --no-daemonize -r "$OMPI_URI_FILE" &
  ompi_server_pid=$!
  sleep 0.5

  if [[ ! -f "$OMPI_URI_FILE" ]]; then
    echo "ERROR: ompi-server did not create URI file" >&2
    kill $ompi_server_pid 2>/dev/null || true
    exit 1
  fi

  MPIRUN_EXTRA_ARGS="--ompi-server file:$OMPI_URI_FILE"
fi

MPIRUN_CMD="$MPIRUN -n 1 $MPIRUN_EXTRA_ARGS"

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM"
MACE_CMD="python3 $MACE_SERVER"

function cleanup {
  kill -9 ${macepid-} ${abinpid-} > /dev/null 2>&1 || true
  if [[ -n "$ompi_server_pid" ]]; then
    kill $ompi_server_pid > /dev/null 2>&1 || true
  fi
  rm -f ompi_uri.txt
  exit 0
}
trap cleanup INT ABRT TERM EXIT

# Launch mock MACE server
$MPIRUN_CMD $MACE_CMD > $MACE_OUT 2>&1 &
macepid=$!

# Wait for the server to write the port file
MAX_WAIT=15
i=0
while [[ ! -f mace_port.txt.1 ]] && [[ $i -lt $MAX_WAIT ]]; do
  sleep 0.5
  let ++i
done

if [[ ! -f mace_port.txt.1 ]]; then
  echo "ERROR: MACE server did not write port file"
  cat $MACE_OUT 2>/dev/null || true
  exit 1
fi

# Launch ABIN
$MPIRUN_CMD $ABIN_CMD > $ABINOUT 2>&1 &
abinpid=$!

# Monitor both processes
MAX_ITER=60
iter=0
while true; do
  abin_alive=$(ps -p $abinpid -o pid= 2>/dev/null | wc -l)
  mace_alive=$(ps -p $macepid -o pid= 2>/dev/null | wc -l)

  if [[ $abin_alive -eq 0 ]] && [[ $mace_alive -eq 0 ]]; then
    break
  elif [[ $abin_alive -eq 0 ]] && [[ $mace_alive -ne 0 ]]; then
    sleep 1
    kill -9 $macepid > /dev/null 2>&1 || true
    break
  elif [[ $mace_alive -eq 0 ]] && [[ $abin_alive -ne 0 ]]; then
    echo "MACE server died. Killing ABIN." >&2
    kill -9 $abinpid > /dev/null 2>&1 || true
    break
  fi

  sleep 0.5
  let ++iter
  if [[ $iter -gt $MAX_ITER ]]; then
    echo "Maximum wait time exceeded."
    break
  fi
done
