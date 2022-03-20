#!/bin/bash

set -euo pipefail

rm -f ERROR ABIN_ERROR? abin.out
if [[ $1 = "clean" ]];then
  exit 0
fi

ABINEXE=$1

# Test that ABIN fails with invalid cmdline argument
$ABINEXE -i input.in -x mini.xyz > abin.out 2>&1 || true
mv ERROR ABIN_ERROR1

# Check that ABIN fails with invalid tau0_langevin
$ABINEXE -i input.in2 -x mini.xyz > abin.out 2>&1 || true
mv ERROR ABIN_ERROR2
