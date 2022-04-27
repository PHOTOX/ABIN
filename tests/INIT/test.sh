#!/bin/bash

set -euo pipefail

rm -f ERROR abin.out movie.xyz restart.xyz *dat
if [[ $1 = "clean" ]];then
  exit 0
fi

ABINEXE=$1

# Testing that ABIN does not run
# if movie.xyz and/or restart.xyz are present
# and irest=0

touch movie.xyz
$ABINEXE -x mini.xyz -i input.in > abin.out 2>&1 || true
mv ERROR ABIN_ERROR1

rm movie.xyz
touch restart.xyz
$ABINEXE -x mini.xyz -i input.in > abin.out 2>&1 || true
mv ERROR ABIN_ERROR2

rm restart.xyz
$ABINEXE -x mini.xyz -i input.in.sh > abin.out 2>&1 || true
mv ERROR ABIN_ERROR3

$ABINEXE -x mini.xyz -i input.in.nhc > abin.out 2>&1 || true
mv ERROR ABIN_ERROR4
