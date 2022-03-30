#/bin/bash

rm -f EXIT ERROR remd.out abin.out restart.xyz.??.old restart.xyz.?? restart.xyz.??.? geom.dat.?? movie.xyz.?? cp2k.out temper.dat.?? energies.dat.??
if [[ "$1" = "clean" ]];then
   exit 0
fi

ABINEXE=$1
ABININ=input.in
ABINVEL=vel0.in
ABINOUT=abin.out

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

N_REPLICAS=$(egrep --only-matching 'nreplica\s*=\s*[0-9]' $ABININ | egrep -o [0-9])

$MPIRUN -np $N_REPLICAS $ABINEXE -i $ABININ -v $ABINVEL > $ABINOUT
$MPIRUN -np $N_REPLICAS $ABINEXE -i ${ABININ}2 >> $ABINOUT

# Test that ABIN stops when file EXIT is present
touch EXIT
$MPIRUN -np $N_REPLICAS $ABINEXE -i ${ABININ}3 >> $ABINOUT
rm -f EXIT

# Useful line in case you need to debug multiple MPI processes
# $MPIRUN -np $N_REPLICAS xterm -e gdb $ABINEXE -i $ABININ -v $ABINVEL > $ABINOUT
