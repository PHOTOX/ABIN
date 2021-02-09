#/bin/bash
set -euo pipefail
# Useful for debugging
#set -x

export ABINEXE=$1

source ../test_tc_server_utils.sh

# The goal here is to test verious failure modes
# and how TC and ABIN responds to them.

# We're going to be testing multiple things
# in this single test, and we will be collecting
# ABIN and TC error messages and compare them with the reference.
# This is not a perfect approach, but let's see how it works.
# Also, Codecov coverage will help us determine
# that we have hit all the paths.

set_default_vars
set_mpich_vars

# If $1 = "clean"; exit early.
if ! clean_output_files $1; then
  exit 0
fi

# Exit early for OpenMPI build.
check_for_openmpi

launch_hydra_nameserver $MPICH_HYDRA

cleanup() {
  kill -9 $hydrapid > /dev/null 2>&1 || true
  exit 0
}

trap cleanup INT ABRT TERM EXIT

# This is used in all individual scripts
# that we call below.
grep_tc_error() {
  tcout=$1
  grep 'what()' $tcout >> $TC_ERROR_FILE
}

export -f grep_tc_error

./test1.sh
./test2.sh

# TODO: Check how ABIN handles MPI error
# (we'll need to build faulty tc_server.

# Check how tc_server handles bad input
# (again, we'll need a modified version)
# Basically we should test every assertion
# in the TCServerMock code.

# Check handling of port.txt file in ABIN.
# (without launching the tc_server)

