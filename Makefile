# Super simple Makefile for ABIN		Daniel Hollas,2014
#
# Simply type "make" and you should get the binary named $OUT
# Before recompiling, it is wise to clean up by "make clean"
#
# The machine-dependent variables are included from make.vars
# No user modification to this file should be necessary.

# WARNING: dependecies on *.mod files are hidden!
# If you change modules, you should recompile the whole thing i.e. make clean;make
include make.vars

export SHELL=/bin/bash
export DATE=`date +"%X %x"`
ifeq ($(shell git --version|cut -b -3),git)
export COMMIT=`git log -1 --pretty=format:"commit %H"`
endif

F_OBJS :=  random.o shake.o nosehoover.o transform.o potentials.o  estimators.o gle.o ekin.o vinit.o  \
nab.o force_bound.o force_guillot.o water.o force_cp2k.o forces.o surfacehop.o force_abin.o  analyze_ext_distp.o density.o analysis.o  \
minimizer.o arrays.o init.o mdstep.o 

C_OBJS := EWALD/ewaldf.o

LIBS = WATERMODELS/libttm.a

ifeq ($(NAB),TRUE)
  C_OBJS += nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o
  LIBS += NAB/libnab.a  NAB/arpack.a  NAB/blas.a
  CFLAGS +=  -INAB/include  
  FFLAGS +=  -DNAB
endif

ifeq ($(FFTW),TRUE)
  LIBS := -lfftw3 ${LIBS}
  FFLAGS := -DUSEFFTW ${FFLAGS}
  F_OBJS := fftw_interface.o ${F_OBJS}
endif

ifeq ($(CP2K),TRUE)
  # The following variables should be the same that were used to compile CP2K.
  # Also , be carefull with FFTW clashes
  FFTWPATH   := /usr/lib/x86_64-linux-gnu/
  FFTW_INC   := /usr/include
  FFTW_LIB   := ${FFTWPATH}
  FFLAGS := -DCP2K -I${FFTW_INC} -I$(BLASPATH)/include ${FFLAGS}
  LDLIBS := -L${CP2KPATH} -lcp2k \
      -L${BLASPATH}/lib $(BLASPATH)/lib/libacml.a \
      ${FFTW_LIB}/libfftw3.a  ${FFTW_LIB}/libfftw3_threads.a\
      ${LDLIBS}
endif

#MPI STUFF
ifeq  ($(MPI),TRUE) 
  FC = $(MPIPATH)/bin/mpif90
  MPIINC = -DMPI -I$(MPIPATH)/include/
  F_OBJS := force_tera.o ${F_OBJS}
endif

LDLIBS = -lm -lstdc++ ${LIBS}

# This hack is needed for force_tera.o and fftw_interface.o
F_OBJS := modules.o utils.o interfaces.o force_mm.o ${F_OBJS}

ALLDEPENDS = ${C_OBJS} ${F_OBJS}


# This is the default target
${OUT} : abin.o
	cd WATERMODELS && make all 
	${FC} ${FFLAGS} WATERMODELS/water_interface.o ${ALLDEPENDS}  $< ${LDLIBS} ${MPILIBS} -o $@

# Always recompile abin.F90 to get current date and commit
abin.o : abin.F90 ${ALLDEPENDS} WATERMODELS/water_interface.cpp
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	$(FC) $(FFLAGS) $(MPIINC) -c abin.F90

clean :
	/bin/rm -f *.o *.mod
	cd WATERMODELS && make clean

distclean :
	/bin/rm -f *.o *.mod NAB/*.o

# Run all tests (this is currently compiler dependent)
# You might expect some machine precision differences
test :
	/bin/bash ./test.sh ${OUT} all
# Test only surface hopping.
testsh :
	/bin/bash ./test.sh ${OUT} sh
# Clean all test folders.
testclean :
	/bin/bash ./test.sh ${OUT} clean

# This will automatically generate new reference data for tests
makeref :
	/bin/bash ./test.sh ${OUT} makeref

# Dummy target for debugging purposes
debug: 
	echo ${LIBS}
	echo ${C_OBJS}
	echo ${CFLAGS}
	echo ${F_OBJS}

.PHONY: clean distclean test testsh testclean makeref debug

.SUFFIXES: .F90 .f90 .f95 .f03 .F03

.F90.o:
	echo "${F_OBJS}"
	$(FC) $(FFLAGS) -c ${MPIINC} $<

.f90.o:
	$(FC) $(FFLAGS) -c ${MPIINC} $<

.f95.o:
	$(FC) $(FFLAGS) -c ${MPIINC} $<

.f03.o:
	$(FC) $(FFLAGS) -c ${MPIINC} $<

.F03.o:
	$(FC) $(FFLAGS) -c ${MPIINC} $<

