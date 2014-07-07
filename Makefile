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
FFLAGS = -g -fopenmp  -Wall -fbounds-check -Og -ffpe-trap=invalid,zero,overflow #static # -O2 -ip -ipo  #-fno-underscoring -fopenmp
CFLAGS =  -g -INAB/include #-Wno-unused-result " 
LD = -lfftw3 -lm -lstdc++


export SHELL=/bin/bash
export DATE=`date +"%X %x"`
export COMMIT=`git log -1 --pretty=format:"commit %H"`

F_OBJS = modules.o interfaces.o random.o nosehoover.o stage.o potentials.o  estimators.o force_mm.o nab.o gle.o analyze_ext_distp.o  \
velverlet.o surfacehop.o minimizer.o force_bound.o respa_shake.o force_guillot.o \
shake.o respa.o analysis.o init.o force_clas.o force_quantum.o density.o ran1.o vinit.o \
shift.o ekin.o force_abin.o

C_OBJS = nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o  EWALD/ewaldf.o

LIBS = NAB/libnab.a  NAB/arpack.a  NAB/blas.a

${OUT} : ${C_OBJS} ${F_OBJS} ${LIBS} date.inc abin.o 
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	${FC} ${FFLAGS} abin.f90 -c    #possibly for the second time,but whatever
	${FC} ${FFLAGS}  ${C_OBJS} ${F_OBJS} ${LIBS} abin.o ${LD} -o $@

#we need this in case date.inc does not exist
date.inc :
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > $@
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> $@

clean :
	/bin/rm -f *.o *.mod NAB/*.o

test :
	/bin/bash ./test.sh ${OUT}

.PHONY: clean test

.SUFFIXES: .F90 .f90 .f95

.F90.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.f95.o:
	$(FC) $(FFLAGS) -c $<

