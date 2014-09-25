#!/bin/bash

if [[ ! -e job.log ]];then
   echo "ERROR: file job.log does not exist. Exiting now..."
   exit 1
fi

NODE=$(head -1 job.log)
JOB=$(tail -1 job.log)

ssh -n $NODE "touch $JOB/EXIT"

