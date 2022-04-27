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
rm -f TC_ERROR? ABIN_ERROR? ${TCOUT}* ${ABINOUT}*
if ! clean_output_files $1; then
  exit 0
fi

# Exit early for OpenMPI build.
check_for_openmpi
# Exit early for IntelMPI.
# IntelMPI does seem to have hydra_nameserver,
# but for some reason it cannot bind to port on GHA.
#check_for_intelmpi

# We relaunch the nameserver in each subtest
# due to the bug in hydra_nameserver
#launch_hydra_nameserver $MPICH_HYDRA

# Compile default TC server
$MPICXX $TCSRC -Wall -o $TCEXE

echo "########### SUBTEST 1 ###################"
./test1.sh
echo "########### SUBTEST 2 ###################"
./test2.sh || true
echo "########### SUBTEST 3 ###################"
./test3.sh
echo "########### SUBTEST 4 ###################"
./test4.sh
# Getting rid of variable MPI trace
head -1 ABIN_ERROR4 > tmp
mv tmp ABIN_ERROR4
echo "########### SUBTEST 5 ###################"
./test5.sh
echo "########### SUBTEST 6 ###################"
./test6.sh
echo "########### SUBTEST 7 ###################"
./test7.sh
echo "########### SUBTEST 8 ###################"
./test8.sh

# TODO: Check how tc_server handles bad input
# (again, we'll need a modified version)
# Basically we should test every assertion
# in the TCServerMock code.
