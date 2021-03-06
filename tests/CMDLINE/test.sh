#!/bin/bash

set -euo pipefail

rm -f ERROR abin.out
if [[ $1 = "clean" ]];then
  exit 0
fi

ABINEXE=$1

# Test that ABIN fails with invalid cmdline argument
$ABINEXE -invalid > abin.out || true

# Test that ABIN prints help, without an error
$ABINEXE -h >> abin.out || echo "ERROR when printing help" >> ERROR
