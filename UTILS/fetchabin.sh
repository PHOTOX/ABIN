#!/bin/bash
# Simple script that fetches data from scratch to the current working directory.
# Expects the presence of file job.log created by r.abin.

if [[ ! -e job.log ]];then
   echo "ERROR: file job.log does not exist. Exiting now..."
   exit 1
fi

NODE=$(head -1 job.log)
JOB=$(tail -1 job.log)

KDE=`pwd`

# copy all data from scratch if it is newer (-u switch)
# and preserve the timestamps (-p switch)
ssh -n $NODE "cp -r -u -p $JOB/* $KDE"

#Try to determine inose form input, ignore errors
inose=$(awk -F"[,=]" '{if($1=="inose")print $2}' $KDE/input.in 2> /dev/null)

#When running microcanonical simulation, check energy conservation

if [[ $inose -eq "0" ]];then
   echo "Checking energy conservation." 
   echo "# Energy_jump   Energy_drift" 
   grep -s -v -e "#" $KDE/energies.dat | checkenergy.awk
fi

