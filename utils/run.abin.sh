#!/bin/bash

# This is a sample script for launching ABIN simulations in cluster environments.
# 1. Copy all files from $PWD to the node's scratch directory.
#    (not copying folders except for the BASH-interface folder)
# 2. Launch ABIN.
# 3. Copy data back (only newer files are copied!).
# 4. Remove scratch directory.

set -euo pipefail

# Example SGE params for PHOTOX clusters
#$ -cwd -notify
#$ -q aq -pe shm 1

# Example SLURM parameters
#SBATCH --mem=1000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

### SETUP ###
JOBNAME=ABIN_${JOB_ID}_$$
ABIN_IN="input.in"
ABIN_OUT="abin.out"
GEOM_IN="mini.xyz"
VELOC_IN=
#############

SCRDIR=/scratch/$USER/$JOBNAME

function copy_to_scrdir() {
  # In SLURM, we might already have a dedicated directory
  if [[ -n ${SLURM_TMPDIR-} ]];then
    SCRDIR=$SLURM_TMPDIR
  elif [[ -d $SCRDIR ]];then
    echo "ERROR: Job directory $SCRDIR already exist!" >&2
    echo "Perhaps it's a leftover from some old job." >&2
    exit 1
  else
    mkdir $SCRDIR
  fi

  # First, copy only files, not directories.
  cp -p $(find . -maxdepth 1 -type f) $SCRDIR/.

  if [[ -d $INTERFACE ]];then
    cp -r $INTERFACE $SCRDIR/
  fi
  if [[ -d $INTERFACE_REF ]];then
    cp -r $INTERFACE_REF $SCRDIR/
  fi
  if [[ -d MM ]];then
    cp -r MM $SCRDIR/
  fi
}

function copy_from_scrdir() {
  echo "Copying from scratch directory back to $LAUNCH_DIR"
  rsync -auvz * $LAUNCH_DIR
  if [[ $? -ne "0" ]];then
    echo "Error when copying the data from scratch directory back to the server." >&2
    echo "I will keep the directory $SCRDIR on node:" >&2
    uname -a >&2
    exit 1
  fi

  rm $LAUNCH_DIR/job.log
  rm -r $SCRDIR
}

function files_exist {
   for file in $* ;
   do
      if [ ! -f $file ];then
         echo "ERROR: Cannot find file $file" >&2
         error=1
      fi
   done
   if [[ -n ${error-} ]];then
      exit 1
   fi
}

function parse_input() {
  # TODO: This is rather brittle, need better regex
  pot=$(awk -F"[! ,=\"']+" '{if($1=="pot")print $2}' $ABIN_IN)
  # Upper case, this is the folder with file interface
  INTERFACE=${pot^^}

  # when using reference potential for multiple timestepping
  pot_ref=$(awk -F"[! ,=\"']+" '{if($1=="pot_ref")print $2}' $ABIN_IN)
  INTERFACE_REF=${pot_ref^^}

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
    REMD=true
    echo "REMD detected. Number of replicas = $N_REPLICAS"
  fi
}

files_exist $ABIN_IN
parse_input

# $ABINEXE exported via SetEnvironment.sh
set +u
source SetEnvironment.sh ABIN
set -u
if [[ ${REMD-} == "true" ]];then
  ABINEXE="$MPIRUN -n $N_REPLICAS $ABINEXE"
fi

uname -n > job.log
echo "$SCRDIR" >> job.log

LAUNCH_DIR=$PWD

# The local SetEnvironment script is sourced in file-based ab initio interfaces
# so that they are system-agnostic.
# We need to use local copy because when global copy
# happened to change, running ABIN simulations were crashing.
# We always copy a new version to local folder,
# since outdated versions can lead to bugs.
# The downside is that when you restart a simulation,
# you are not quite guaranteed that you will use the same version
# of the ab initio program (e.g. when new version is installed on your cluster).
OURSETENV=$(which SetEnvironment.sh)
cp $OURSETENV .

copy_to_scrdir

trap copy_from_scrdir EXIT SIGUSR2

cd $SCRDIR

if [[ -z ${VELOC_IN-} ]];then
   $ABINEXE -i $ABIN_IN -x $GEOM_IN >> $ABIN_OUT
else
   $ABINEXE -v $VELOC_IN -i $ABIN_IN -x $GEOM_IN >> $ABIN_OUT
fi
