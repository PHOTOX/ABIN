#!/bin/bash
NODE=$(head -1 job.log)
JOB=$(tail -1 job.log)

ssh -n $NODE "touch $JOB/EXIT"

