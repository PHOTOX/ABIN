# Super simple Makefile for ABIN		Daniel Hollas,2014
#
# Simply type "make" and you should get the binary named $OUT
# Before recompiling, it is wise to clean up by "make clean"

# WARNING: dependecies on *.mod files are hidden!
# If you change modules, you should recompile the whole thing i.e. make clean;make

OUT = abin.dev
#FC = /usr/local/programs/common/intel/compiler/2013.5.192/bin/ifort
FC = gfortran
CC = gcc
# if you have FFTW libraries available, set it to TRUE
# if not, some ABIN functionality will be limited
FFTW = TRUE 
# Should we compile with NAB libraries (AMBER force field)
# Currently only possible with gfortran
NAB  = TRUE

# -----------------------------------------------------------------------

FFLAGS :=  -g  -fopenmp # -Wall -Wextra -fbounds-check -ffpe-trap=invalid,zero,overflow #static # -O2 -ip -ipo  #-fno-underscoring -fopenmp
CFLAGS :=   -g #-Wno-unused-result " 

export SHELL=/bin/bash
export DATE=`date +"%X %x"`
ifeq ($(shell git --version|cut -b -3),git)
export COMMIT=`git log -1 --pretty=format:"commit %H"`
endif

F_OBJS := utils.o interfaces.o random.o shake.o nosehoover.o transform.o potentials.o  estimators.o gle.o ekin.o vinit.o  \
force_mm.o nab.o force_bound.o force_guillot.o water.o forces.o surfacehop.o force_abin.o  analyze_ext_distp.o density.o analysis.o  \
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

LDLIBS = -lm -lstdc++ ${LIBS}

F_OBJS := modules.o ${F_OBJS}

ALLDEPENDS = ${C_OBJS} ${F_OBJS}


# This is the default target
${OUT} : abin.o
	cd WATERMODELS && make all 
	${FC} ${FFLAGS} WATERMODELS/water_interface.o ${ALLDEPENDS}  $< ${LDLIBS} -o $@

# Always recompile abin.F90 to get current date and commit
abin.o : abin.F90 ${ALLDEPENDS} WATERMODELS/water_interface.cpp
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	$(FC) $(FFLAGS) -c abin.F90

clean :
	/bin/rm -f *.o *.mod

cleanall :
	/bin/rm -f *.o *.mod NAB/*.o
	cd WATERMODELS && make clean

# Run all tests (this is currently compiler dependent)
# You might expect some machine precision differences
test :
	/bin/bash ./test.sh ${OUT} all
# Test only surface hopping.
testsh :
	/bin/bash ./test.sh ${OUT} sh
# Clean all test folders.
testcl :
	/bin/bash ./test.sh ${OUT} clean

# This will automatically generate new reference data for tests
makeref :
	/bin/bash ./test.sh ${OUT} makeref

# Dummy target for debugging purposes
debug: 
	echo ${LIBS}
	echo ${C_OBJS}
	echo ${CFLAGS}

.PHONY: clean test testsh testcl makeref debug

.SUFFIXES: .F90 .f90 .f95 .f03 .F03

.F90.o:
	echo "${F_OBJS}"
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.f95.o:
	$(FC) $(FFLAGS) -c $<

.f03.o:
	$(FC) $(FFLAGS) -c $<

.F03.o:
	$(FC) $(FFLAGS) -c $<

