#!/bin/bash

# Install OpenMPI implementation of Message Passing Interface (MPI)
# This script installs both the MPI compilers (mpifort,mpicc)
# and MPI process manager (mpirun)

# Exit script immediately upon error
set -euo pipefail

# Path as an optional first parameter
PLUMED_DIR="${1-/home/$USER/plumed}"
# We take the current stable version as default
# (as of 17 Dec 2020)
PLUMED_VERSION="${2-"2.6.2"}"
# Use two processors for the build
NPROC=2

TAR_FILE="plumed-src-${PLUMED_VERSION}.tgz"
DOWNLOAD_URL="https://github.com/plumed/plumed2/releases/download/v${PLUMED_VERSION}/${TAR_FILE}"
INSTALL_DIR="$PLUMED_DIR/$PLUMED_VERSION/install"

if [[ -d $PLUMED_DIR/$PLUMED_VERSION ]];then
  echo "Found existing Plumed installation in $PLUMED_DIR/$PLUMED_VERSION"
  echo "Remove this folder if you want to reinstall"
  exit 1
fi

mkdir -p $PLUMED_DIR/$PLUMED_VERSION/src
mkdir -p $PLUMED_DIR/$PLUMED_VERSION/pkg

curl -L "$DOWNLOAD_URL" > $PLUMED_DIR/$PLUMED_VERSION/pkg/${TAR_FILE}
cd $PLUMED_DIR/$PLUMED_VERSION/src && tar -xzf ../pkg/${TAR_FILE} 
cd plumed-${PLUMED_VERSION}

# To keep things simple, we don't compile with MPI support
# We also disable the use of external BLAS and LAPACK
# to prevent Gfortran version conflicts 
# According to the manual, Plumed should ship with their default 
# BLAS and LAPACK, that will be fine for us for now.

# --enable-static-patch is needed by ABIN integration
# It is the default, but let's be explicit
./configure CXX=g++ CC=gcc \
  --disable-mpi --disable-xdrfile \
  --disable-external-lapack --disable-external-blas \
  --enable-static-patch --enable-shared \
  --prefix=${INSTALL_DIR} 2>&1 |\
  tee configure.log

make -j $NPROC 2>&1 | tee make.log

make install 2>&1 | tee make_install.log

echo "
Succesfully installed PLUMED!
Set the following in your ABIN make.vars

PLUMED = TRUE
PLUMED_INC = ${INSTALL_DIR}/lib/plumed/src/lib/Plumed.inc

or rerun configure as

./configure --plumed ${INSTALL_DIR}
"
