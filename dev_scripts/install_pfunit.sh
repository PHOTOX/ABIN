#!/bin/bash

# Exit script immediately upon error
set -euo pipefail

REPO_DIR="$HOME/pfunit"
if [[ "$#" -eq 1 && ! -z $1 ]];then
   REPO_DIR=$1
fi

if [[ -e $REPO_DIR ]];then
  echo "ERROR: $REPO_DIR already exists."
  exit 1
fi

git clone --recursive https://github.com/Goddard-Fortran-Ecosystem/pFUnit $REPO_DIR && cd $REPO_DIR

# This seems to be the last commit that works with GCC-7
# Later commits fail with internal compiler error when compiling fArgParse submodule. 
git checkout --recurse-submodules 0a09db354b665f1518e36460396c348c19185e04

mkdir -p build && cd build
export FC=gfortran
export FFLAGS=

# TODO: Make it possible to change final install dir via shell parameter
cmake .. -DGIT_SUBMODULE=OFF -DSKIP_FHAMCREST=YES
make tests
make install

echo "
Succesfully installed pFUnit library 
Set the following path in your make.vars

PFUNIT_PATH = $REPO_DIR/build/installed/

or rerun configure as
./configure --pfunit $REPO_DIR/build/installed/
"
