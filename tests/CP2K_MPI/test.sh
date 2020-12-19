#/bin/bash

if [[ $1 = "clean" ]];then
   rm -f restart.xyz.old output *.diff temper.dat energies.dat restart.xyz WATER* movie.xyz cp2k.out
   exit 0
fi

ABINEXE=$1

# $MPI_PATH exported at the top of  test.sh
$MPI_PATH/bin/mpirun  -np 3 $ABINEXE -i input.in > output
