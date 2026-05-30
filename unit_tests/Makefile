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
  FFLAGS = $(PFUNIT_EXTRA_FFLAGS) -Wall -Wno-maybe-uninitialized -O0 -ffpe-trap=invalid,zero,overflow,denormal
  FFLAGS += -I../src/
  FFLAGS += $(DFLAGS)
endif

utils_TESTS := test_utils.pf
utils_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
utils_OTHER_SRCS := throw_with_pfunit.F90
utils_EXTRA_USE := throw_with_pfunit_mod
utils_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,utils))

utils_driver.o: throw_with_pfunit.o

masses_TESTS := test_masses.pf
masses_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
masses_OTHER_SRCS := throw_with_pfunit.F90
masses_EXTRA_USE := throw_with_pfunit_mod
masses_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,masses))

prng_TESTS := test_prng.pf
prng_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
prng_OTHER_SRCS := throw_with_pfunit.F90
prng_EXTRA_USE := throw_with_pfunit_mod
prng_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,prng))

# This is needed for compatibility with new Plumed 2.8
# See: https://github.com/plumed/plumed2/pull/757
export PLUMED_DP2CUTOFF_NOSTRETCH = 1

plumed_TESTS := test_plumed.pf
plumed_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
plumed_OTHER_SRCS := throw_with_pfunit.F90
plumed_EXTRA_USE := throw_with_pfunit_mod
plumed_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,plumed))

terapi_TESTS := test_terapi.pf
terapi_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
terapi_OTHER_SRCS := throw_with_pfunit.F90
terapi_EXTRA_USE := throw_with_pfunit_mod
terapi_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,terapi))

mace_TESTS := test_mace.pf
mace_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
mace_OTHER_SRCS := throw_with_pfunit.F90
mace_EXTRA_USE := throw_with_pfunit_mod
mace_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,mace))

fftw_TESTS := test_fftw.pf
fftw_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
fftw_OTHER_SRCS := throw_with_pfunit.F90
fftw_EXTRA_USE := throw_with_pfunit_mod
fftw_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,fftw))

gle_TESTS := test_gle.pf
gle_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
gle_OTHER_SRCS := throw_with_pfunit.F90
gle_EXTRA_USE := throw_with_pfunit_mod
gle_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,gle))

nhc_TESTS := test_nhc.pf
nhc_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
nhc_OTHER_SRCS := throw_with_pfunit.F90
nhc_EXTRA_USE := throw_with_pfunit_mod
nhc_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,nhc))

pot_TESTS := test_potentials.pf
pot_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
pot_OTHER_SRCS := throw_with_pfunit.F90
pot_EXTRA_USE := throw_with_pfunit_mod
pot_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,pot))

mm_TESTS := test_mm.pf
mm_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
mm_OTHER_SRCS := throw_with_pfunit.F90
mm_EXTRA_USE := throw_with_pfunit_mod
mm_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,mm))

io_TESTS := test_io.pf
io_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
io_OTHER_SRCS := throw_with_pfunit.F90
io_EXTRA_USE := throw_with_pfunit_mod
io_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,io))

spline_TESTS := test_spline.pf
spline_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
spline_OTHER_SRCS := throw_with_pfunit.F90
spline_EXTRA_USE := throw_with_pfunit_mod
spline_EXTRA_INITIALIZE := initialize_throw
# TODO: Re-enable this when splines are re-implemented
#$(eval $(call make_pfunit_test,spline))

interface_TESTS := test_interface.pf
interface_OTHER_LIBRARIES := ${LDFLAGS} ${LDLIBS}
interface_OTHER_SRCS := throw_with_pfunit.F90
interface_EXTRA_USE := throw_with_pfunit_mod
interface_EXTRA_INITIALIZE := initialize_throw
$(eval $(call make_pfunit_test,interface))

# TODO: Re-enable spline tests
#all : utils prng masses plumed fftw terapi gle nhc pot mm spline interface io
all : utils prng masses plumed fftw terapi mace gle nhc pot mm interface io

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

test : clean_plumed_output clean_terapi_output clean_gle_output
	./utils
	./prng
	./masses
	./plumed
	./fftw
	./terapi
	./mace
	./gle
	./nhc
	./pot
	./mm
	#./spline
	./io
	./interface

clean : clean_plumed_output clean_terapi_output clean_gle_output
	/bin/rm -f *dat *.o *.out *.mod *.inc test_*.F90 ERROR
	/bin/rm -f masses utils plumed fftw terapi mace gle pot mm spline interface

clean_plumed_output :
	/bin/rm -f bck.* plumed_*.dat

clean_terapi_output :
	/bin/rm -f wfn.bin*

clean_gle_output :
	/bin/rm -f GLE-S.bin GLE-T.bin
 
.PHONY: clean test clean_plumed_output clean_terapi_output

.SUFFIXES: .F90
