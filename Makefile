#!/bin/bash
#setting variables needed for icc
export LANG=C
export LC_ALL=C 

NAB=true #true or false
#FFLAGS=" -O0 -ffpe-trap=invalid,zero,overflow -g -static"   #DEBUG
#FFLAGS="-fopenmp -O2 -pg"
#CFLAGS="-pg -O2 -pthread"  #PARALLEL VERSION
FFLAGS="  -g  "  # -ffpe-trap=invalid,zero,overflow -g  "  #-O2 -ip -ipo " #-fno-underscoring -fopenmp"
CFLAGS=" -g " #-Wno-unused-result " 
OUT=abin.dev
FCC=gfortran
CC=gcc

rm *.o ../BUILD/$OUT

$FCC $FFLAGS  -c modules.f90  nosehoover.f90 stage.f90 estimators.f90  nab.F90 gle.F90 analyze_ext_distp.f90 potentials.f90 velverlet.f90 surfacehop.f90 force_mm.f90 minimizer.f90 random.f force_bound.f90  respa_shake.f90 force_guillot.f90 shake.f90 abin.f90 respa.f90 analysis.f90 init.F90 force_clas.f90 force_quantum.f90 density.f90 ran1.f vinit.f shift.f90 ekin.f90 force_abin.f90

if [ "$NAB" = "true" ] ;then
$CC -c $CFLAGS -INAB/include nabinit_pme.c NAB/sff_my_pme.c NAB/memutil.c NAB/prm.c NAB/nblist_pme.c NAB/binpos.c  EWALD/ewaldf.c 
#$FCC $FFLAGS -c nab.F90
fi

$FCC $FFLAGS  *.o NAB/libnab.a  NAB/arpack.a  NAB/blas.a -lfftw3 -lm -lstdc++ -o BUILD/$OUT

rm *.o 

