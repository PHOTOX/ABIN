# Super simple Makefile for ABIN		Daniel Hollas,2014
#
# Simply type "make" and you should get the binary named $OUT
# Before recompiling, it is wise to clean up by "make clean"

# WARNING: dependecies on *.mod files are hidden!
# If you change modules, you should recompile the whole thing i.e. make clean;make

OUT = abin.dev
# You actually have to use gfortran and gcc, because of precompiled NAB libraries
FC = gfortran
CC = gcc
# if you have FFTW libraries available, set it to TRUE
# if not, some ABIN functionality will be limited
FFTW =TRUE 
CP2K =TRUE
BLASPATH = -L/usr/local/lib/acml5.3.1/gfortran64/lib -lacml
CP2KPATH = /usr/local/src/cp2k-2.6.1/lib/Linux-x86-64-gfortran/sopt/

# -----------------------------------------------------------------------

FFLAGS := -g
#FFLAGS :=  -g -fopenmp  -Wall -Wextra -fbounds-check -ffpe-trap=invalid,zero,overflow #static # -O2 -ip -ipo  #-fno-underscoring -fopenmp
CFLAGS :=  -g -INAB/include #-Wno-unused-result " 

LIBS = NAB/libnab.a  NAB/arpack.a  NAB/blas.a WATERMODELS/libttm.a
LDLIBS = -lm -lstdc++ ${LIBS}



export SHELL=/bin/bash
export DATE=`date +"%X %x"`
ifeq ($(shell git --version|cut -b -3),git)
export COMMIT=`git log -1 --pretty=format:"commit %H"`
endif

F_OBJS := utils.o interfaces.o random.o shake.o nosehoover.o transform.o potentials.o  estimators.o gle.o ekin.o vinit.o  \
force_mm.o nab.o force_bound.o force_guillot.o water.o force_cp2k.o forces.o surfacehop.o force_abin.o  analyze_ext_distp.o density.o analysis.o  \
minimizer.o arrays.o init.o mdstep.o 

C_OBJS := nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o  EWALD/ewaldf.o

ifeq ($(FFTW),TRUE)
LDLIBS := -lfftw3 ${LDLIBS}
FFLAGS := -DUSEFFTW ${FFLAGS}
F_OBJS := fftw_interface.o ${F_OBJS}
endif

ifeq ($(CP2K),TRUE)
LDLIBS := -L${CP2KPATH} ${BLASPATH} -lcp2k -lcp2kbase  -lfftw3 ${LDLIBS}
#LDLIBS := ${LDLIBS} -lcp2ksubsys -lcp2kstart -lcp2kcommon -lcp2kinput -lcp2kmain -lcp2kmachine -lcp2kmotion -lcp2kma -lcp2kao -lcp2kfm -lcp2kgrid \
	-lcp2kmpiwrap -lcp2kpw -lcp2ktmc -lcp2kxc
FFLAGS := -DCP2K ${FFLAGS}
#F_OBJS := ${F_OBJS} force_cp2k.o 
endif

F_OBJS := modules.o ${F_OBJS}

ALLDEPENDS = ${C_OBJS} ${F_OBJS}

# This is the default target
${OUT} : abin.o
	cd WATERMODELS && make all 
	${FC} ${FFLAGS} WATERMODELS/water_interface.o ${ALLDEPENDS}  $< ${LDLIBS} -o $@

#Always recompile abin.F03 to get current date and commit
abin.o : abin.F03 ${ALLDEPENDS} WATERMODELS/water_interface.cpp
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	$(FC) $(FFLAGS) -c abin.F03

clean :
	/bin/rm -f *.o *.mod

cleanall :
	/bin/rm -f *.o *.mod NAB/*.o
	cd WATERMODELS && make clean

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

