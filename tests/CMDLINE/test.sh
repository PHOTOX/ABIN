#!/bin/bash

set -euo pipefail

rm -f ERROR abin.out
if [[ $1 = "clean" ]];then
  exit 0
fi

ABINEXE=$1

# Test that ABIN fails with invalid cmdline argument
$ABINEXE -invalid > abin.out 2>&1 || true
mv ERROR ABIN_ERROR1

# Test that ABIN prints help, without an error
$ABINEXE -h >> abin.out 2>&1 || echo "ERROR when printing help" >> ABIN_ERROR1
$ABINEXE --help >> abin.out 2>&1 || echo "ERROR when printing help" >> ABIN_ERROR1

# Test that ABIN prints version, without an error
$ABINEXE -V >> abin.out 2>&1 || echo "ERROR when printing version" >> ABIN_ERROR1
$ABINEXE --version >> abin.out 2>&1 || echo "ERROR when printing version" >> ABIN_ERROR1

# Test that ABIN does not accept empty arguments
$ABINEXE -v >> abin.out 2>&1 || true
mv ERROR ABIN_ERROR2

$ABINEXE -x >> abin.out 2>&1 || true
mv ERROR ABIN_ERROR3

$ABINEXE -i >> abin.out 2>&1 || true
mv ERROR ABIN_ERROR4

$ABINEXE -M >> abin.out 2>&1 || true
mv ERROR ABIN_ERROR5
