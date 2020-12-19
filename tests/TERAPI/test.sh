#/bin/bash
#module load mpich2/1.4.1p1
#MPIRUN=mpiexec

MPIRUN=$MPI_PATH/bin/mpirun

rm -f restart.xyz movie.xyz
if [[ "$1" = "clean" ]];then
   rm -f  output terapi.out temper.dat energies.dat
   exit 0
fi



#DHHack:
if [[ -z $1 ]];then
   ABINEXE=../../abin
else
   ABINEXE=$1
fi

$MPIRUN  -np 1 ./tera-mpiapi > terapi.out &
# Get PID of the last process
terapid=$!

sleep 2
# Ugly workaround because MPI_Lookup does not work
grep port_name: terapi.out | awk '{print $6}' > port.txt
$MPIRUN  -np 1 $ABINEXE > output &
abinpid=$!

while true;do
   sleep 1
   if ! `ps|grep -q $terapid` && ! `ps|grep -q $abinpid` ;then
      echo "Both ABIN and TeraChem stopped."
      break
   fi
   if ! `ps|grep -q $terapid` ;then
      echo "Terachem died. Killing ABIN."
      kill -9 $abinpid 
      break
   fi   
   if ! `ps|grep -q $abinpid` ;then
      echo "ABIN died. Killing TeraChem."
      echo $terapid
      kill -9 $terapid 
      break
   fi
done

