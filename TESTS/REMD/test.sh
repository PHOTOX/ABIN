#/bin/bash
source SetEnvironment.sh CP2K

rm -f  movie.xyz* cp2k.out temper.dat.* energies.dat.* fort.*
rm -f restart.xyz*

$MPIRUN  -np 2 ../../abin.cp2k.popt
#$MPIRUN -np 2 xterm -e gdb ../../abin.mpi
