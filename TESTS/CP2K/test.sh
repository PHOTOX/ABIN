#/bin/bash
source SetEnvironment.sh CP2K

rm -f restart.xyz WATER* movie.xyz cp2k.out

$MPIRUN -np 4 ../../abin.cp2k.popt > output
