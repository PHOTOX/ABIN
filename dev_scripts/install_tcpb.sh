#!/bin/bash
set -euo pipefail

REPO_URL="https://github.com/mtzgroup/tcpb-cpp.git"
REPO_DIR="$HOME/tcpb-cpp"
if [[ "$#" -eq 1 && ! -z $1 ]];then
   REPO_DIR=$1
fi

if [[ -e $REPO_DIR ]];then
  echo "ERROR: $REPO_DIR already exists."
  exit 1
fi

git clone --recursive "${REPO_URL}" "${REPO_DIR}" && cd "$REPO_DIR"

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
