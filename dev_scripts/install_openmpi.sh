#!/bin/bash

# Install OpenMPI implementation of Message Passing Interface (MPI)
# This script installs both the MPI compilers (mpifort)
# and MPI process manager (mpirun)

# Exit script immediately upon error
set -euo pipefail

# Path as an optional first parameter
OPENMPI_DIR="${1-/home/$USER/openmpi}"
# We take current stable version as default
# (as of 06 Nov 2020)
OPENMPI_VERSION=${2-"4.0"}
OPENMPI_VERSION_PATCH=${3-"0"}

# TODO: Make this variable
# Github Actions machines have two CPUs, per:
# https://docs.github.com/en/free-pro-team@latest/actions/reference/specifications-for-github-hosted-runners#supported-runners-and-hardware-resources
NCPUS=2

TAR_FILE="openmpi-${OPENMPI_VERSION}.${OPENMPI_VERSION_PATCH}.tar.gz"
DOWNLOAD_URL="https://download.open-mpi.org/release/open-mpi/v${OPENMPI_VERSION}/${TAR_FILE}"
INSTALL_DIR="$OPENMPI_DIR/$OPENMPI_VERSION/install"

if [[ -d "$OPENMPI_DIR/$OPENMPI_VERSION" ]];then
  echo "Found existing OPENMPI installation in $OPENMPI_DIR/$OPENMPI_VERSION"
  echo "Remove this folder if you want to reinstall"
  exit 1
fi

mkdir -p $OPENMPI_DIR/$OPENMPI_VERSION/src
mkdir -p $OPENMPI_DIR/$OPENMPI_VERSION/pkg

curl "${DOWNLOAD_URL}" > $OPENMPI_DIR/$OPENMPI_VERSION/pkg/${TAR_FILE}
cd $OPENMPI_DIR/$OPENMPI_VERSION/src && tar -xzf ../pkg/${TAR_FILE} && \
  cd openmpi-${OPENMPI_VERSION}.${OPENMPI_VERSION_PATCH}

# For compilation on PHOTOX clusters, add --with-sge
./configure FC=gfortran CC=gcc CXX=g++ --with-gnu-ld --prefix=${INSTALL_DIR} 2>&1 |\
  tee configure.log
make -j ${NCPUS} all 2>&1 | tee make.log
make install 2>&1 | tee make_install.log

echo "
Succesfully installed OpenMPI!
Set the following path in your ABIN make.vars

MPI_PATH = ${INSTALL_DIR}

or rerun configure as

./configure --mpi ${INSTALL_DIR}
"
