#/bin/bash

# The goal here is to test various failure modes
# and how TC and ABIN respond to them.

# We're testing multiple things in this single test,
# and we will be collecting ABIN and TC error messages
# and compare them with the reference.
# This is not a perfect approach, but let's see how it works.
# Also, Codecov coverage will help us determine
# that we have hit all the paths.
set -euo pipefail

export ABINEXE=$1

source ../test_tc_server_utils.sh

set_default_vars
set_mpich_vars

# If $1 = "clean"; exit early.
rm -f TC_ERROR? ${TCOUT}* ${ABINOUT}*
if ! clean_output_files $1; then
  exit 0
fi

# Exit early for OpenMPI build.
check_for_openmpi

launch_hydra_nameserver $MPICH_HYDRA

# Compile default TC server
$MPICXX $TCSRC -Wall -o $TCEXE

cleanup() {
  #kill -9 $hydrapid > /dev/null 2>&1 || true
  exit 0
}

trap cleanup INT ABRT TERM EXIT

./test1.sh
cat abin.out1 tc.out1

./test2.sh
./test3.sh
echo "########### TEST 3 ###################"
cat abin.out3 tc.out3.1

# TODO: Check how ABIN handles MPI error
# (we'll need to build faulty tc_server.

# Check how tc_server handles bad input
# (again, we'll need a modified version)
# Basically we should test every assertion
# in the TCServerMock code.

# Check handling of port.txt file in ABIN.
# (without launching the tc_server)
