#!/bin/bash
# This is a sample script for launching an ABIN simulation
# in cluster environments using MPI TeraChem interface.

set -euo pipefail

# SGE Params on PHOTOX clusters
#$ -V -cwd -notify

# ABIN SETUP 
ABIN_OUT=abin.out
ABIN_IN=input.in
GEOM_IN=mini.xyz
VELOC_IN=

TC_IN=tc.inp
MPITYPE=0   # 0 - ground state AIMD
            # 2 - Surface Hopping
################

JOBNAME=TERABIN_${JOB_ID}_$$

USE_HYDRA_NAMESERVER="false"

SCRDIR=/scratch/$USER/$JOBNAME
LAUNCH_DIR=$PWD

function copy_to_scrdir() {
  if [[ -d $SCRDIR ]];then
    echo "ERROR: Job directory $SCRDIR already exist!" >&2
    echo "Perhaps it's a leftover from some old job." >&2
    exit 1
  fi

  mkdir $SCRDIR

  # Copying only files, not directories.
  cp -p $(find . -maxdepth 1 -type f) $SCRDIR/.

  if [[ -d $INTERFACE_REF ]];then
    cp -r $INTERFACE_REF $SCRDIR/
  fi
}

function copy_from_scrdir() {
  echo "Copying scratch back to $LAUNCH_DIR"
  rsync -auvz * $LAUNCH_DIR
  if [[ $? -ne "0" ]];then
    echo "Error when copying the data from scratch back to the server." >&2
    echo "I will keep the directory $SCRDIR on node:" >&2
    uname -a >&2
    exit 1
  fi

  rm $LAUNCH_DIR/job.log
  rm -r $SCRDIR
}

function files_exist() {
   for file in $* ;
   do
      if [[ ! -f $file ]];then
         echo "ERROR: Cannot find file $file" >&2
         error=1
      fi
   done
   if [[ -n ${error-} ]];then
      exit 1
   fi
}

function launch_hydra_nameserver() {
  hydra=$(ps -C hydra_nameserver -o pid= || true)
  if [[ -z $hydra ]]; then
    MPIPATH=$(dirname $MPIRUN)
    echo "Launching hydra nameserver for MPI_Lookup"
    if [[ -f $MPIPATH/hydra_nameserver ]]; then
      $MPIPATH/hydra_nameserver &
    else
      echo "ERROR: Could not find hydra_nameserver executable" >&2
      exit 1
    fi
  fi
}

function grep_tc_port {
  tcout=$1
  portfile=$2
  maxiter=10
  i=0
  # Grep port name from TC output, pass to ABIN via a file.
  tcport=$(awk '/port_name:/ {print $(NF);exit}' $tcout)
  while [[ -z ${tcport} ]]; do
    if [[ $i -gt $maxiter ]];then
      echo "ERROR: Could not grep port name from $tcout" >&2
      exit 1
    fi
    sleep 1
    tcport=$(awk '/port_name:/ {print $(NF);exit}' $tcout)
    let ++i
  done
  echo $tcport > $portfile
}

function validate_inputs() {
  files_exist $TC_IN $ABIN_IN

  # Check pot=_tera_ in ABIN input
  test=$(egrep -o -e ^[^!]*pot[[:space:]]*=[[:space:]]*[\'\"]_tera_[\"\'] $ABIN_IN || true)
  if [[ -z $test ]];then
    echo "ERROR: You did not specify pot='_tera_' in $ABIN_IN." >&2
    exit 1 
  fi
}

function parse_inputs() {

  N_TERA_SERVERS=$(egrep -o -e ^[^!]*nteraservers[[:space:]]*=[[:space:]]*[0-9]+\\b $ABIN_IN | egrep -o -e [0-9]+|| true)
  # Set default if not specified in the ABIN input.
  : ${N_TERA_SERVERS:=1}

  # automatic REMD settings from ABIN input
  N_REPLICAS=1
  test=$(egrep -o -e ^[^!]*iremd[[:space:]]*=[[:space:]]*1\\b $ABIN_IN || true)
  if [[ -n $test ]];then
    echo "REMD detected. Assuming one TC server per replica!"
    N_REPLICAS=$(egrep -o -e ^[^!]*nreplica[[:space:]]*=[[:space:]]*[0-9]+\\b $ABIN_IN | egrep -o -e [0-9]+|| true)
    if [[ -z $N_REPLICAS ]];then
      echo "ERROR: Could not determine the number of replicas from $ABIN_IN" >&2
      exit 1
    fi
    N_TERA_SERVERS=$N_REPLICAS
  fi

  # When using a reference potential in multiple time-step MD,
  # we need to copy an extra directory.
  pot_ref=$(awk -F"[! ,=\"']+" '{if($1=="pot_ref")print $2}' $ABIN_IN)
  INTERFACE_REF=${pot_ref^^}

  if [[ -n ${NSLOTS-} && $NSLOTS -ne $N_TERA_SERVERS ]];then
    echo "ERROR: Number of requested GPUs ($NSLOTS)" >&2
    echo "is not equal to the number of TC servers ($N_TERA_SERVERS)" >&2
    exit 1
  fi
}

# All auxiliary functions defined, let's GO!
validate_inputs
parse_inputs

# SET THE ENVIRONMENT
set +u
source SetEnvironment.sh TERACHEM
source SetEnvironment.sh ABIN
set -u

if [[ $N_TERA_SERVERS -gt 1 ]];then
  # hydra_nameserver <= v3.4.1 does not work with multiple servers.
  USE_HYDRA_NAMESERVER=false
fi

export OMP_NUM_THREADS=$N_TERA_SERVERS
MPIRUN_ABIN="$MPIRUN -n $N_REPLICAS"
MPIRUN_TERA="$MPIRUN -n 1"
if [[ $USE_HYDRA_NAMESERVER = "true" ]];then
  MPIRUN_ABIN="${MPIRUN_ABIN} -nameserver $HOSTNAME"
  MPIRUN_TERA="${MPIRUN_TERA} -nameserver $HOSTNAME"
fi

# HERE WE GO!
uname -n > job.log
echo "$SCRDIR" >> job.log

copy_to_scrdir
cd $SCRDIR

TC_PORT=$JOBNAME.$$
if [[ $USE_HYDRA_NAMESERVER = "true" ]]; then
  launch_hydra_nameserver
fi

trap copy_from_scrdir EXIT SIGUSR2

declare -A job_pids
# LAUNCH TERACHEM
for ((itera=1;itera<=N_TERA_SERVERS;itera++)) {
   gpuid=$((itera-1))
   $MPIRUN_TERA $TERAEXE -g$gpuid --inputfile=$TC_IN --UseMPI=$MPITYPE --MPIPort=$TC_PORT.$itera > $TC_IN.out.$itera 2>&1 &
   job_pids[tc_$itera]=$!
}

if [[ $USE_HYDRA_NAMESERVER != "true" ]];then
  TC_PORT_FILE="port.txt"
  # Grep port names from TC output, pass to ABIN via a file.
  for ((itera=1;itera<=N_TERA_SERVERS;itera++)); do
    grep_tc_port $TC_IN.out.$itera  $TC_PORT_FILE.$itera
  done
fi

# LAUNCH ABIN
ABIN_CMD="$ABINEXE -i $ABIN_IN -x $GEOM_IN"
if [[ $USE_HYDRA_NAMESERVER = "true" ]];then
  ABIN_CMD="$ABIN_CMD -M $TC_PORT"
fi
if [[ -n $VELOC_IN ]];then
   ABIN_CMD=$ABIN_CMD" -v $VELOC_IN"
fi
$MPIRUN_ABIN $ABIN_CMD > $ABIN_OUT 2>&1 &
job_pids[abin]=$!

# PERIODICALLY CHECK WHETHER ABIN AND TC ARE RUNNING
function join_by { local IFS="$1"; shift; echo "$*"; }
regex=`join_by \| ${job_pids[@]}`
while true;do
   njobs=$(ps -eo pid | egrep "$regex" | wc -l)
   if [[ $njobs -eq 0 ]];then
      echo "Both ABIN and TeraChem stopped"
      break
   elif [[ $njobs -lt ${#job_pids[@]} ]];then
      echo "One of the TC servers or ABIN died. Killing the rest." >&2
      kill -9 ${job_pids[@]} 
      break
   fi
   sleep 3
done
