#!/bin/bash

# Install OpenMPI implementation of Message Passing Interface (MPI)
# This script installs both the MPI compilers (mpifort,mpicc)
# and MPI process manager (mpirun)

# Exit script immediately upon error
set -eou pipefail

# Path as an optional first parameter
MPICH_DIR="${1-/home/$USER/mpich}"
# We take current stable version as default
# (as of 06 Nov 2020)
MPICH_VERSION="${2-"3.3.2"}"

TAR_FILE="mpich-${MPICH_VERSION}.tar.gz"
DOWNLOAD_URL="https://www.mpich.org/static/downloads/${MPICH_VERSION}/${TAR_FILE}"
INSTALL_DIR="$MPICH_DIR/$MPICH_VERSION/install"

if [[ -d $MPICH_DIR/$MPICH_VERSION ]];then
  echo "Found existing MPICH installation in $MPICH_DIR/$MPICH_VERSION"
  echo "Remove this folder if you want to reinstall"
  exit 1
fi

mkdir -p $MPICH_DIR/$MPICH_VERSION/src
mkdir -p $MPICH_DIR/$MPICH_VERSION/pkg

curl "$DOWNLOAD_URL" > $MPICH_DIR/$MPICH_VERSION/pkg/${TAR_FILE}
cd $MPICH_DIR/$MPICH_VERSION/src && tar -xzf ../pkg/${TAR_FILE} && cd mpich-${MPICH_VERSION}

# If you're building MPI for general use, not only for ABIN,
# you might want to change the configure options:
# --enable-fortran=yes Compile all versions of Fortran interfaces
#                     This option is needed for GitHub Actions build, I don't know why
# --with-hydra-pm=pmiserv --with-namepublisher=pmi
#                     Needed for MPI interface with TeraChem
./configure FC=gfortran CC=gcc \
  --enable-fortran=all --disable-cxx \
  --with-hydra-pm=pmiserv --with-namepublisher=pmi \
  --enable-static --disable-shared \
  --prefix=${INSTALL_DIR} 2>&1 |\
  tee configure.log
make 2>&1 | tee make.log
make install 2>&1 | tee make_install.log

echo "
Succesfully installed MPICH!
Set the following path in your ABIN make.vars

MPI_PATH = ${INSTALL_DIR}

or rerun configure as

./configure --mpi ${INSTALL_DIR}
"
