FFLAGS= -fopenmp -g #-O2  #PARALLEL VERSION
#CFLAGS="-pg -O2 -pthread"  #PARALLEL VERSION
#FFLAGS =  -g -Wall -fbounds-check #-O0 -ffpe-trap=invalid,zero,overflow -g static "  #-O2 -ip -ipo " #-fno-underscoring -fopenmp"
CFLAGS =  -g -INAB/include #-Wno-unused-result " 
OUT = abin.dev.openmp
FC = gfortran
CC = gcc
LD = -lfftw3 -lm -lstdc++


F_OBJS = modules.o nosehoover.o stage.o estimators.o nab.o gle.o analyze_ext_distp.o potentials.o \
velverlet.o surfacehop.o force_mm.o minimizer.o random.f force_bound.o respa_shake.o force_guillot.o \
shake.o abin.o respa.o analysis.o init.o force_clas.o force_quantum.o density.o ran1.o vinit.o \
shift.o ekin.o force_abin.o

C_OBJS = nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o  EWALD/ewaldf.o

LIBS = NAB/libnab.a  NAB/arpack.a  NAB/blas.a

abin.dev : ${C_OBJS} ${F_OBJS} ${LIBS}
	${FC} ${FFLAGS}  ${C_OBJS} ${F_OBJS} ${LIBS} ${LD} -o $@

clean :
	/bin/rm *.o *.mod NAB/*.o

.SUFFIXES: .F90 .f90

.F90.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

