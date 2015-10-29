#/bin/bash
source SetEnvironment.sh CP2K

rm -f restart.xyz WATER* movie.xyz cp2k.out

$MPIRUN -np 3 ../../abin.remd
#$MPIRUN -np 2 xterm -e gdb ../../abin.cp2k.popt
