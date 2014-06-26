# Makefile for ABIN		Daniel Hollas,2014
# Simply type "make" and you should get the binary
#
# WARNING:dependecies on *.mod files are hidden!
# if you change modules, you should recompile the whole thing i.e. make clean;make
#
OUT = abin.dev
FC = gfortran
CC = gcc

#CFLAGS="-pg -O2 -pthread"  #PARALLEL VERSION
FFLAGS =  -g -Wall -fbounds-check -O0 -ffpe-trap=invalid,zero,overflow  #static # -O2 -ip -ipo " #-fno-underscoring -fopenmp"
CFLAGS =  -g -INAB/include #-Wno-unused-result " 
LD = -lfftw3 -lm -lstdc++


export SHELL=/bin/bash
export DATE=`date +"%X %x"`
export COMMIT=`git log -1 --pretty=format:"commit %H"`

F_OBJS = modules.o nosehoover.o stage.o estimators.o nab.o gle.o analyze_ext_distp.o potentials.o \
velverlet.o surfacehop.o force_mm.o minimizer.o random.f force_bound.o respa_shake.o force_guillot.o \
shake.o abin.o respa.o analysis.o init.o force_clas.o force_quantum.o density.o ran1.o vinit.o \
shift.o ekin.o force_abin.o

C_OBJS = nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o  EWALD/ewaldf.o

LIBS = NAB/libnab.a  NAB/arpack.a  NAB/blas.a

${OUT} : ${C_OBJS} ${F_OBJS} ${LIBS}
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	${FC} ${FFLAGS} abin.f90 -c    #possibly for the second time,but whatever
	${FC} ${FFLAGS}  ${C_OBJS} ${F_OBJS} ${LIBS} ${LD} -o $@

clean :
	/bin/rm -f *.o *.mod NAB/*.o

.PHONY: clean

.SUFFIXES: .F90 .f90 .f95

.F90.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.f95.o:
	$(FC) $(FFLAGS) -c $<

