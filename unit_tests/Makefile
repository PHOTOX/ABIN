# Makefile for ABIN Unit Tests
# Should be invoked from top level Makefile,
# which defines certain variables used here,
# namely path to PFUNIT library
 
# Compile Unit Tests using pFUnit library
#
# TODO: Exit prematurely if PFUNIT_PATH is not defined
ifneq ($(strip $(PFUNIT_PATH)),)
  LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_PATH)/PFUNIT-4.*) | xargs -n1 | sort -V))
  include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk
  # TODO: Here we rewrite FFLAGS, but we should be reusing it
  # Need to figure out how to do test coverage separately
  FFLAGS = $(PFUNIT_EXTRA_FFLAGS)
  FFLAGS += -I../
endif
 
STATIC_LIBS = ../libabin.a ../WATERMODELS/libwater.a

LDLIBS = -lm -lstdc++ ${LIBS}

utils_TESTS := test_utils.pf
utils_REGISTRY :=
utils_OTHER_SOURCES :=
utils_OTHER_LIBRARIES := -lgcov -fprofile-arcs -L"../" -labin -L"../WATERMODELS/" -lwater $(LDLIBS)
utils_OTHER_INCS :=

$(eval $(call make_pfunit_test,utils))

all : utils

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

test :
	./utils

clean :
	/bin/rm -f *.o *.mod *.inc test_*.F90 utils
 
.PHONY: clean test

.SUFFIXES: .F90

.F90.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -c $<

