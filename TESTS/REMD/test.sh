#/bin/bash

rm -f remd.out restart.xyz.??.old output output2 restart.xyz.?? restart.xyz.??.? geom.dat.?? movie.xyz.?? cp2k.out temper.dat.?? energies.dat.?? 
if [[ "$1" = "clean" ]];then
   exit 0
fi

ABINEXE=$1
MPIRUN="$MPI_PATH/bin/mpirun"
# If this fails, just try plain mpirun and hope for the best
if [[ ! -f $MPIRUN ]];then
   MPIRUN=mpirun
fi

$MPIRUN  -np 3 $ABINEXE -i input.in -v vel0.in > output
$MPIRUN  -np 3 $ABINEXE -i input.in2 >> output
#$MPIRUN -np 2 xterm -e gdb ../../abin.mpi
