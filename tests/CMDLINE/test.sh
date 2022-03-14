#!/bin/bash

set -euo pipefail

rm -f ERROR abin.out
if [[ $1 = "clean" ]];then
  exit 0
fi

ABINEXE=$1

# Test that ABIN fails with invalid cmdline argument
$ABINEXE -invalid > abin.out 2>&1 || true

# Test that ABIN prints help, without an error
$ABINEXE -h >> abin.out 2>&1 || echo "ERROR when printing help" >> ERROR
$ABINEXE --help >> abin.out 2>&1 || echo "ERROR when printing help" >> ERROR

# Test that ABIN prints version, without an error
$ABINEXE -V >> abin.out 2>&1 || echo "ERROR when printing version" >> ERROR
$ABINEXE --version >> abin.out 2>&1 || echo "ERROR when printing version" >> ERROR
