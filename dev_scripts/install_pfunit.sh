#!/bin/bash

# Exit script immediately upon error
set -euo pipefail

# This must correspond to a git tag in https://github.com/Goddard-Fortran-Ecosystem/pFUnit/tags
PFUNIT_VERSION="v4.18.0"

if [[ "$#" -ne 1 ]];then
  echo "Provide path where to install pfunit as a first parameter"
  exit 1
fi

mkdir -p "$1"
REPO_DIR="$1"

git clone --recursive --branch $PFUNIT_VERSION https://github.com/Goddard-Fortran-Ecosystem/pFUnit "$REPO_DIR" && cd "$REPO_DIR"

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
