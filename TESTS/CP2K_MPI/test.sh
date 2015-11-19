#/bin/bash
source SetEnvironment.sh CP2K

if [[ $1 = "clean" ]];then
   rm -f restart.xyz.old output *.diff temper.dat energies.dat restart.xyz WATER* movie.xyz cp2k.out
   exit 0
fi

ABINEXE=$1

$MPIRUN  -np 3 $ABINEXE -i input.in > output
