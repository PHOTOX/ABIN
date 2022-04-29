# Makefile for ABIN

# The user-defined variables are included from file "make.vars',
# which is not under version control
# No user modification to this Makefile file should be necessary.

# Type "make" and you should get the binary named src/$BIN as defined in make.vars
# Before recompiling, it is wise to clean up by "make clean"

# WARNING: dependecies on *.mod files are not properly resolved here!
# If you change modules, you should recompile the whole thing by running
# $ make clean && make

# By default run all end-to-end tests in tests/
TEST=all
# ABIN binary name
BIN=abin
# Optional compilation parameters
# Some functionality will not work without them
MPI=FALSE
FFTW=FALSE
CP2K=FALSE
PLUMED=FALSE
# https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html
LDLIBS=
LDFLAGS=

# Export all vars into submake commands
export
# User-defined compilation parameters are in make.vars
# and should override the defaults defined above
include make.vars

export SHELL=/bin/bash
export COMMIT=NaN
ifeq ($(shell git --version | cut -b -3),git)
  ifneq ($(shell git status 2>&1 | head -1 | cut -b -5),fatal)
    export COMMIT=`git rev-parse HEAD`
  endif
endif

ifeq ($(strip $(FFTW)),TRUE)
  DFLAGS += -DUSE_FFTW
  ifeq ($(CP2K),TRUE)
    $(info "!!!!!-------------WARNING---------------!!!!!!!")
    $(info "Using FFTW flag with CP2K may lead to troubles!")
    $(info "!!!!!-------------WARNING---------------!!!!!!!")
    $(info "")
  else
    LDLIBS += -lfftw3
  endif
endif

ifeq ($(strip $(CP2K)),TRUE)
  DFLAGS += -DUSE_CP2K
  # OpenMP does not work with -fno-underscoring
  FFLAGS += -fno-underscoring -fno-openmp
  # The following variables should be the same that were used to compile CP2K.
  # Also, be carefull with FFTW clashes
  LDLIBS += -lcp2k ${CP2K_LIBS}
  LDFLAGS += -L${CP2K_PATH}
endif

ifeq ($(strip $(PLUMED)),TRUE)
 include ${PLUMED_INC}
 DFLAGS += -DUSE_PLUMED
 LDLIBS += ${PLUMED_STATIC_LOAD}
endif

ifeq  ($(strip $(MPI)),TRUE) 
  DFLAGS += -DUSE_MPI
endif

LDLIBS = -labin -lwater $(LDLIBS) -lm -lstdc++
LDFLAGS += -fopenmp -L../src/ -L../water_potentials/

# This is the default target
${BIN} :
	mkdir -p bin/
	$(MAKE) -C water_potentials all
	$(MAKE) -C src $(BIN)
	$(MAKE) -C utils all

utils: ${BIN}
	$(MAKE) -C utils all

clean:
	/bin/rm -f *.gcov bin/*
	$(MAKE) -C water_potentials clean
	$(MAKE) -C src clean
	$(MAKE) -C utils clean
ifneq ($(strip $(PFUNIT_PATH)),)
	$(MAKE) -C unit_tests clean
endif

# Build and run Unit tests (powered by pFUnit library)
unittest: ${BIN}
ifneq ($(strip $(PFUNIT_PATH)),)
	$(MAKE) -C unit_tests all
	$(MAKE) -C unit_tests test
else
	echo "pFUnit library not available. Skipping unit tests."
endif

# End-To-End (E2E) tests
e2etest: ${BIN}
	/bin/bash tests/test.sh ${BIN} $(TEST) ${MPI} ${FFTW} ${PLUMED} ${CP2K} test

# Runs both end-to-end and unit tests
test: unittest e2etest

# Clean all test folders.
testclean:
ifneq ($(strip $(PFUNIT_PATH)),)
	$(MAKE) -C unit_tests clean
endif
	/bin/bash tests/test.sh ${BIN} $(TEST) ${MPI} ${FFTW} $(PLUMED) ${CP2K} clean

# This will automatically generate new reference data for E2E tests
makeref: ${BIN}
	/bin/bash tests/test.sh ${BIN} $(TEST) ${MPI} ${FFTW} $(PLUMED) ${CP2K} makeref


.PHONY: clean test testclean makeref unittest e2etest
