#!/bin/bash
# This is a sample script for launching an ABIN simulation
# in cluster environments using Protocol Buffers (PB) TeraChem interface.

set -uo pipefail

# SGE Params on PHOTOX clusters
#$ -V -cwd -notify

# ABIN SETUP 
ABIN_OUT=abin.out
ABIN_IN=input.in
GEOM_IN=mini.xyz
VELOC_IN=

TC_IN=tc.inp
TCPB_CPP_LIB= # path_to/tcpb-cpp/lib/
################
if [[ -n $TCPB_CPP_LIB ]];then
  export LD_LIBRARY_PATH=$TCPB_CPP_LIB:${LD_LIBRARY_PATH-}
fi

JOBNAME=TCPBABIN_${JOB_ID-}_$$

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
  kill $abin_pid $tcpb_pid &> /dev/null
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

function validate_inputs() {
  files_exist $TC_IN $ABIN_IN

  # Check pot=_tera_ in ABIN input
  test=$(egrep -o -e ^[^!]*pot[[:space:]]*=[[:space:]]*[\'\"]_tcpb_[\"\'] $ABIN_IN || true)
  if [[ -z $test ]];then
    echo "ERROR: You did not specify pot='_tcpb_' in $ABIN_IN." >&2
    exit 1 
  fi
}

function parse_inputs() {
  # automatic REMD settings from ABIN input
  N_REPLICAS=1
  test=$(egrep -o -e ^[^!]*iremd[[:space:]]*=[[:space:]]*1\\b $ABIN_IN || true)
  if [[ -n $test ]];then
    echo "ERROR: REMD is currently not supported with TCPB interface."
    exit 1
  fi

  # When using a reference potential in multiple time-step MD,
  # we need to copy an extra directory.
  pot_ref=$(awk -F"[! ,=\"']+" '{if($1=="pot_ref")print $2}' $ABIN_IN)
  INTERFACE_REF=${pot_ref^^}
}

# All auxiliary functions defined, let's GO!
validate_inputs
parse_inputs

# SET THE ENVIRONMENT
set +u
source SetEnvironment.sh TERACHEM
source SetEnvironment.sh ABIN
set -u

uname -n > job.log
echo "$SCRDIR" >> job.log

copy_to_scrdir
cd $SCRDIR

trap copy_from_scrdir EXIT SIGUSR2

# $RANDOM returns 0-32767, which is within port range on Linux.
# We add 1025 to avoid clashing with system-reserved ports.
TC_PORT=$((1025+RANDOM))
TC_HOST=localhost
TC_OUT=$(basename $TC_IN .inp)".tcpb.out"

#### LAUNCH TERACHEM ####
$TERAEXE -s $TC_PORT &> $TC_OUT &
tcpb_pid=$!
# 10 seconds is a conservative number for TCPB startup
sleep 10

ABIN_CMD="$ABINEXE -i $ABIN_IN -x $GEOM_IN"
if [[ -n $VELOC_IN ]];then
   ABIN_CMD=$ABIN_CMD" -v $VELOC_IN"
fi
ABIN_CMD=$ABIN_CMD" --tcpb-host $TC_HOST --tcpb-port $TC_PORT --tcpb-input-file $TC_IN"

#### LAUNCH ABIN ####
$ABIN_CMD &> $ABIN_OUT &
abin_pid=$!

# PERIODICALLY CHECK WHETHER ABIN AND TC ARE RUNNING
while true;do
   njobs=$(ps -eo pid | egrep "$tcpb_pid|$abin_pid" | wc -l)
   if [[ $njobs -eq 0 ]];then
      echo "Both ABIN and TeraChem stopped"
      break
   elif [[ $njobs -lt 2 ]];then
      echo "TCPB server or ABIN died. Killing the other." >&2
      kill $abin_pid $tcpb_pid &> /dev/null
      break
   fi
   sleep 5
done
