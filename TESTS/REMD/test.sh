#/bin/bash
source SetEnvironment.sh CP2K

if [[ "$1" = "clean" ]];then
   rm -f remd.out restart.xyz.??.old output output2 restart.xyz.?? restart.xyz.??.? geom.dat.?? movie.xyz.?? cp2k.out temper.dat.?? energies.dat.?? 
   exit 0
fi

ABINEXE=$1

$MPIRUN  -np 3 $ABINEXE -i input.in -v vel0.in > output
$MPIRUN  -np 3 $ABINEXE -i input.in2 >> output
#$MPIRUN -np 2 xterm -e gdb ../../abin.mpi
