# Makefile for ABIN		Daniel Hollas,2014
# Simply type "make" and you should get the binary
# Before recompiling, it is wise to clean up by "make clean"
#
# WARNING:dependecies on *.mod files are hidden!
# if you change modules, you should recompile the whole thing i.e. make clean;make
#
OUT = abin.dev
FC = gfortran
CC = gcc

#CFLAGS="-pg -O2 -pthread"  #PARALLEL VERSION
FFLAGS =  -g  -Wall -Wextra -fbounds-check -Og -ffpe-trap=invalid,zero,overflow #static # -O2 -ip -ipo  #-fno-underscoring -fopenmp
CFLAGS =  -g -INAB/include #-Wno-unused-result " 
LIBS = NAB/libnab.a  NAB/arpack.a  NAB/blas.a
LDLIBS = -lfftw3 -lm -lstdc++ ${LIBS}


export SHELL=/bin/bash
export DATE=`date +"%X %x"`
export COMMIT=`git log -1 --pretty=format:"commit %H"`

F_OBJS = modules.o utils.o interfaces.o random.o shake.o nosehoover.o stage.o potentials.o  estimators.o gle.o ekin.o vinit.o  \
force_mm.o nab.o analyze_ext_distp.o velverlet.o surfacehop.o minimizer.o force_bound.o respa_shake.o force_guillot.o \
respa.o density.o analysis.o init.o force_clas.o force_quantum.o  \
shift.o force_abin.o

C_OBJS = nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o  EWALD/ewaldf.o

ALLDEPENDS = ${C_OBJS} ${F_OBJS}

${OUT} : abin.o
	${FC} ${FFLAGS} ${ALLDEPENDS} $< ${LDLIBS} -o $@

#Always recompile abin.f90 to get current date and commit
abin.o : abin.f90 ${ALLDEPENDS}
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	$(FC) $(FFLAGS) -c abin.f90

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

.SUFFIXES: .F90 .f90 .f95

.F90.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.f95.o:
	$(FC) $(FFLAGS) -c $<

