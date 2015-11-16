# Super simple Makefile for ABIN		Daniel Hollas,2014
#
# The user defined variables are included from file make.vars,
# which is not under version control
#
# Simply type "make" and you should get the binary named $OUT
# Before recompiling, it is wise to clean up by "make clean"
#
# The machine-dependent variables are included from make.vars
# No user modification to this file should be necessary.

# WARNING: dependecies on *.mod files are hidden!
# If you change modules, you should recompile the whole thing i.e. make clean;make
#
# For compilation with static system libraries, see:
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539
#
TEST=all
include make.vars

export SHELL=/bin/bash
export DATE=`date +"%X %x"`
ifeq ($(shell git --version|cut -b -3),git)
export COMMIT=`git log -1 --pretty=format:"commit %H"`
endif

F_OBJS :=  transform.o potentials.o  estimators.o gle.o ekin.o vinit.o  \
force_nab.o force_bound.o force_guillot.o water.o force_cp2k.o forces.o surfacehop.o force_abin.o  analyze_ext_distp.o density.o analysis.o  \
minimizer.o arrays.o init.o mdstep.o 

C_OBJS := EWALD/ewaldf.o

LIBS += WATERMODELS/libttm.a

ifeq ($(NAB),TRUE)
  C_OBJS += nabinit_pme.o NAB/sff_my_pme.o NAB/memutil.o NAB/prm.o NAB/nblist_pme.o NAB/binpos.o
  LIBS   += NAB/libnab.a  NAB/arpack.a  NAB/blas.a
  CFLAGS +=  -INAB/include  
  DFLAGS +=  -DNAB
endif

ifeq ($(FFTW),TRUE)
  ifneq ($(CP2K),TRUE)
   LIBS := -lfftw3 ${LIBS}
  endif
  DFLAGS += -DUSEFFTW
  F_OBJS := fftw_interface.o ${F_OBJS}
endif

ifeq ($(CP2K),TRUE)
  # The following variables should be the same that were used to compile CP2K.
  # Also, be carefull with FFTW clashes
  DFLAGS += -DCP2K 
  LIBS += -L${CP2KPATH} -lcp2k ${CP2K_LIBS} 
endif

#MPI STUFF
ifeq  ($(MPI),TRUE) 
  DFLAGS += -DMPI
  INC    += $(MPI_INC)
  LIBS   += $(MPI_LIBS)
  F_OBJS := force_tera.o remd.o ${F_OBJS}
endif

LDLIBS = -lm -lstdc++ ${LIBS}
# The following line does not seem to work
#LDLIBS = ${LIBS} -static-libgfortran -Wl,-Bstatic -lstdc++ -lm -Wl,-Bdynamic  

# Adding rest of the Fortran objects
# This hack is needed for force_tera.o and fftw_interface.o
F_OBJS := modules.o utils.o interfaces.o force_mm.o random.o shake.o nosehoover.o  ${F_OBJS}

ALLDEPENDS = ${C_OBJS} ${F_OBJS}


# This is the default target
${OUT} : abin.o
	cd WATERMODELS && make all 
	${FC} ${FFLAGS} WATERMODELS/water_interface.o ${ALLDEPENDS} $< ${LDLIBS} -o $@

# Always recompile abin.F90 to get current date and commit
abin.o : abin.F90 ${ALLDEPENDS} WATERMODELS/water_interface.cpp
	echo "CHARACTER (LEN=*), PARAMETER :: date ='${DATE}'" > date.inc
	echo "CHARACTER (LEN=*), PARAMETER :: commit='${COMMIT}'" >> date.inc
	$(FC) $(FFLAGS) $(DFLAGS) $(INC) -c abin.F90

clean :
	/bin/rm -f *.o *.mod $(OUT)
	cd WATERMODELS && make clean

# Remove NAB objects as well
distclean :
	/bin/rm -f *.o *.mod NAB/*.o

# Run test suite 
test :
	/bin/bash TESTS/test.sh ${OUT} $(TEST) ${NAB} ${MPI} ${CP2K} ${FFTW}

# Clean all test folders.
testclean :
	/bin/bash TESTS/test.sh ${OUT} clean 

# This will automatically generate new reference data for tests
makeref :
	/bin/bash TESTS/test.sh ${OUT} $(TEST) ${NAB} ${MPI} ${CP2K} ${FFTW} makeref
 
# Dummy target for debugging purposes
debug: 
	echo ${LIBS}
	echo ${INC}
	echo ${DFLAGS}
	echo ${CFLAGS}
	echo ${FFLAGS}
	echo ${C_OBJS}
	echo ${F_OBJS}

.PHONY: clean distclean test testsh testclean makeref debug

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

