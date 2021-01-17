#/bin/bash

rm -f remd.out restart.xyz.??.old output output2 restart.xyz.?? restart.xyz.??.? geom.dat.?? movie.xyz.?? cp2k.out temper.dat.?? energies.dat.?? 
if [[ "$1" = "clean" ]];then
   exit 0
fi

ABINEXE=$1

# If MPI_PATH is not set, let's hope mpirun is in PATH
if [[ ! -d $MPI_PATH ]];then
   MPIRUN=mpirun
elif [[ -f "$MPI_PATH/bin/orterun" ]];then
  # OpenMPI does not allow oversubscribing by default,
  # i.e. more MPI processes than CPU cores.
  # Github Actions runners have only 2 CPUs,
  # so we need to explicitly set oversubscription
  MPIRUN="$MPI_PATH/bin/orterun -oversubscribe"
else
  MPIRUN="$MPI_PATH/bin/mpirun"
fi


$MPIRUN -np 3 $ABINEXE -i input.in -v vel0.in > abin.out
$MPIRUN -np 3 $ABINEXE -i input.in2 >> abin.out
#$MPIRUN -np 2 xterm -e gdb ../../abin.mpi
