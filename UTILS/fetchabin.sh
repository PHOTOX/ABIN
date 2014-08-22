#!/bin/bash
# Simple script that fetches data from scratch to the current working directory.
# Expects the presence of file job.log created by r.abin.

NODE=$(head -1 job.log)
JOB=$(tail -1 job.log)

KDE=`pwd`

ssh -n $NODE "cp -p $JOB/* $KDE"

#Try to determine inose form input, ignore errors
inose=$(awk -F"[,=]" '{if($1=="inose")print $2}' $KDE/input.in 2> /dev/null)
if [[ $inose -eq "0" ]];then
   echo "Checking energy conservation." 
   echo "# Energy_jump   Energy_drift" 
   grep -s -v -e "#" $KDE/energies.dat | checkenergy.awk
fi

