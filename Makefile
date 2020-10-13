# Super simple Makefile for ABIN		Daniel Hollas,2014

# The user defined variables are included from file "make.vars',
# which is not under version control
# No user modification to this Makefile file should be necessary.

# Simply type "make" and you should get the binary named $BIN
# Before recompiling, it is wise to clean up by "make clean"

# WARNING: dependecies on *.mod files are hidden!
# If you change modules, you should recompile the whole thing i.e. make clean;make

# For compilation with static system libraries, see:
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539

TEST=all
BIN=abin
include make.vars

export SHELL=/bin/bash
export DATE=`date +"%X %x"`

# TODO: Make this if check less stupid, we need to actually 
# check we're in Git repo
ifeq ($(shell git --version|cut -b -3),git)
export COMMIT=`git log -1 --pretty=format:"commit %H"`
endif


F_OBJS := arrays.o transform.o potentials.o estimators.o gle.o ekin.o vinit.o plumed.o \
          force_bound.o water.o force_cp2k.o sh_integ.o surfacehop.o landau_zener.o\
          force_tera.o force_terash.o force_abin.o en_restraint.o analyze_ext_template.o density.o analysis.o \
          minimizer.o mdstep.o forces.o abin.o


LIBS += WATERMODELS/libttm.a

ifeq ($(strip $(FFTW)),TRUE)
  ifneq ($(CP2K),TRUE)
   LIBS := -lfftw3 ${LIBS}
  endif
  DFLAGS += -DUSEFFTW
  F_OBJS := fftw_interface.o ${F_OBJS}
endif

ifeq ($(strip $(CP2K)),TRUE)
  DFLAGS += -DCP2K 
  FFLAGS += -fno-underscoring 
  # OpenMP does not work with -fno-underscoring
  FFLAGS += -fno-openmp
  # The following variables should be the same that were used to compile CP2K.
  # Also, be carefull with FFTW clashes
  LIBS += -L${CP2K_PATH} -lcp2k ${CP2K_LIBS} 
ifeq ($(strip $(FFTW)),TRUE)
   $(info "!!!!!-------------WARNING---------------!!!!!!!")
   $(info "Using FFTW flag with CP2K may lead to troubles!")
   $(info "!!!!!-------------WARNING---------------!!!!!!!")
   $(info "")
endif
endif

ifeq ($(strip $(PLUMED)),TRUE)
 include ${PLUMED_LINK}
 DFLAGS += -DPLUM
 LIBS += ${PLUMED_STATIC_LOAD}
endif

ifeq  ($(strip $(MPI)),TRUE) 
  DFLAGS += -DMPI
  INC    += $(MPI_INC)
  LIBS   += $(MPI_LIBS)
  F_OBJS := remd.o ${F_OBJS}
endif

# Compile Unit Tests using pFUnit library
ifeq  ($(strip $(PFUNIT)),TRUE)
  LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_PATH)/PFUNIT-4.*) | xargs -n1 | sort -V))
  include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk
  FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
endif


LDLIBS = -lm -lstdc++ ${LIBS}
# The following line does not seem to work
#LDLIBS = ${LIBS} -static-libgfortran -Wl,-Bstatic -lstdc++ -lm -Wl,-Bdynamic  

# Adding rest of the Fortran objects
# This hack is needed for force_tera.o and fftw_interface.o
F_OBJS := modules.o utils.o fortran_interfaces.o io.o force_mm.o random.o shake.o nosehoover.o  ${F_OBJS}

# DH: Not sure why this ifeq is neccessary
ifeq ($(strip $(CP2K)),TRUE)
ALLDEPENDS = ${F_OBJS}
else
ALLDEPENDS = ${F_OBJS} WATERMODELS/water_interface.o
endif


# This is the default target
${BIN} : init.o
	cd WATERMODELS && make all 
	${FC} ${FFLAGS} ${ALLDEPENDS} $< ${LDLIBS} -o $@

# Always recompile init.F90 to get current date and commit
# TODO: Figure out a cleaner way to do this
init.o : init.F90 ${ALLDEPENDS}
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -DDATE="'${DATE}'" -DCOMMIT="'${COMMIT}'" -c init.F90

clean :
	cd WATERMODELS && make clean
	/bin/rm -f *.o *.mod libabin.a $(BIN)

# Run the test suite
# TODO: Pass MPI_PATH as well
test : ${BIN}
	/bin/bash TESTS/test.sh ${BIN} $(TEST) ${MPI} ${FFTW} ${PLUM} ${CP2K}

# Clean all test folders.
testclean :
	/bin/bash TESTS/test.sh ${BIN} clean ${MPI} ${FFTW} $(PLUM) ${CP2K}

# This will automatically generate new reference data for tests
makeref :
	/bin/bash TESTS/test.sh ${BIN} $(TEST) ${MPI} ${FFTW} $(PLUM) ${CP2K} makeref

libabin.a : init.o
	ar cru libabin.a init.o $(ALLDEPENDS) && ranlib libabin.a

unittest : libabin
	# TODO: compile and execute unittests

 
# Dummy target for debugging purposes
debug: 
	echo ${LIBS}
	echo ${INC}
	echo ${DFLAGS}
	echo ${CFLAGS}
	echo ${FFLAGS}
	echo ${F_OBJS}

.PHONY: clean distclean test testclean makeref debug

.SUFFIXES: .F90 .f90 .f95 .f03 .F03

.F90.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -c $<

.f90.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -c $<

.f95.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -c $<

.f03.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -c $<

.F03.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -c $<

