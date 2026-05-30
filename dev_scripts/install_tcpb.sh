#!/bin/bash
set -euo pipefail

REPO_URL="https://github.com/mtzgroup/tcpb-cpp.git"
if [[ "$#" -ne 1 ]];then
  echo "Provide path where to install pfunit as a first parameter"
  exit 1
fi
REPO_DIR=$1

if [[ -e $REPO_DIR ]];then
  echo "ERROR: $REPO_DIR already exists."
  exit 1
fi

git clone "${REPO_URL}" "${REPO_DIR}" && cd "$REPO_DIR"

export FC=gfortran
export FFLAGS=

./configure --prefix="$REPO_DIR" gnu
make install

echo "
Succesfully installed TCPB-CPP library 
Set the following path in your make.vars

TCPB_LIB = $REPO_DIR/lib

or rerun configure as
./configure --tcpb $REPO_DIR/lib
"
