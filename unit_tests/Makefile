# Makefile for ABIN Unit Tests
#
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
  FFLAGS = $(PFUNIT_EXTRA_FFLAGS) -Wall -O0 -ffpe-trap=invalid,zero,overflow,denormal
  FFLAGS += -I../src/
  FFLAGS += $(DFLAGS)
endif

# NOTE: $(LIBS) are defined in the root Makefile

utils_TESTS := test_utils.pf
utils_OTHER_LIBRARIES := -lgcov -fprofile-arcs -L"../src/" -labin -L"../water_potentials/" -lwater ${LIBS}
utils_OTHER_SRCS := throw_with_pfunit.F90
utils_EXTRA_USE := throw_with_pfunit_mod
utils_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,utils))

utils_driver.o: throw_with_pfunit.o

plumed_TESTS := test_plumed.pf
plumed_OTHER_LIBRARIES := -lgcov -fprofile-arcs -L"../src/" -labin -L"../water_potentials/" -lwater ${LIBS}
$(eval $(call make_pfunit_test,plumed))

fftw_TESTS := test_fftw.pf
fftw_OTHER_LIBRARIES := -lgcov -fprofile-arcs -L"../src/" -labin -L"../water_potentials/" -lwater ${LIBS}
$(eval $(call make_pfunit_test,fftw))

all : utils plumed fftw

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

test : clean_plumed_output
	./utils
	./plumed
	./fftw

clean : clean_plumed_output
	/bin/rm -f *dat *.o *.out *.mod *.inc test_*.F90 utils plumed fftw

clean_plumed_output :
	/bin/rm -f bck.* plumed_*.dat
 
.PHONY: clean test

.SUFFIXES: .F90
