# Super simple Makefile for ABIN		Daniel Hollas,2014
# Simply type "make" and you should get the binary named $OUT
# Before recompiling, it is wise to clean up by "make clean"
#
# WARNING: dependecies on *.mod files are hidden!
# If you change modules, you should recompile the whole thing i.e. make clean;make
#
OUT = abin.dev
# You actually have to use gfortran and gcc, becouse of precompiled LIBS
FC = gfortran
CC = gcc

FFLAGS =  -g -fopenmp  -Wall -Wextra -fbounds-check -Og -ffpe-trap=invalid,zero,overflow #static # -O2 -ip -ipo  #-fno-underscoring -fopenmp
CFLAGS =  -g -INAB/include #-Wno-unused-result " 
#CFLAGS="-pg -O2 -pthread"  #PARALLEL VERSION

LIBS = NAB/libnab.a  NAB/arpack.a  NAB/blas.a
LDLIBS = -lfftw3 -lm -lstdc++ ${LIBS}

export SHELL=/bin/bash
export DATE=`date +"%X %x"`
export COMMIT=`git log -1 --pretty=format:"commit %H"`

F_OBJS = modules.o utils.o interfaces.o random.o shake.o nosehoover.o stage.o potentials.o  estimators.o gle.o ekin.o vinit.o  \
force_mm.o nab.o force_bound.o force_guillot.o  forces.o surfacehop.o force_abin.f90  analyze_ext_distp.o density.o analysis.o  \
minimizer.o arrays.o init.o mdstep.o 

C_OBJS = nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o  EWALD/ewaldf.o

ALLDEPENDS = ${C_OBJS} ${F_OBJS}

${OUT} : abin.o
	${FC} ${FFLAGS} ${ALLDEPENDS} $< ${LDLIBS} -o $@

#Always recompile abin.f90 to get current date and commit
abin.o : abin.f03 ${ALLDEPENDS}
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	$(FC) $(FFLAGS) -c abin.f03

clean :
	/bin/rm -f *.o *.mod NAB/*.o

test :
	/bin/bash ./test.sh ${OUT} all
testsh :
	/bin/bash ./test.sh ${OUT} sh
testcl :
	/bin/bash ./test.sh ${OUT} clean

makeref :
	/bin/bash ./test.sh ${OUT} makeref

.PHONY: clean test testsh testcl makeref

.SUFFIXES: .F90 .f90 .f95 .f03 .F03

.F90.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.f95.o:
	$(FC) $(FFLAGS) -c $<

.f03.o:
	$(FC) $(FFLAGS) -c $<

.F03.o:
	$(FC) $(FFLAGS) -c $<

