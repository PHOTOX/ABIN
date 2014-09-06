#!/bin/bash

#Main script for Wigner initial conditions
#See README for more detail

#---SETUP--------------------------------
i=1   	    #number of first geometry
ntrajs=1    #total number of generated  geometries
natom=10    #number of atoms
writemass=0 #=1 for atoms other than H, C, N, O, S, F, Cl: see below
iseed=10066077 #random seed 
#-----------------------------------------

inputdeck=input.com

MyIRandom $iseed $ntrajs > iran.dat
if [[ $? -ne "0" ]];then
   echo "Error during random number generation.Exiting now..."
   exit 1
fi


#!DO NOT! EDIT THE FOLLOWING ######
mkdir -p FMSTRAJS
m12root=$(readlink -f /usr/local/programs/molpro/molpro2012.1/arch/amd64-intel_12.0.5.220/molpros_2012_1_Linux_x86_64_i8)
MOLPROEXE=$m12root/bin/molpro

while [ $i -le $ntrajs ]
do
   echo ' &control'> FMSINPOUT/Control.dat
   head -$i iran.dat|tail -1|awk -F" " '{print "IRndSeed=",$1}' >> FMSINPOUT/Control.dat
   echo "
 NumParticles=$natom
 InitState=1
 numstates=1
 IMethod=0	! AIMD, no FMS
 InitialCond='Wigner'
 TimeStep=10.0
 SimulationTime=1.0  
 zEDatFile=.false.   !write Energy file?
 zNDatFile=.false.   !write Population file?
 zTrajFile=.true.   !write Trajectory file?
 /">> FMSINPOUT/Control.dat

# EDIT THE FOLLOWING IF NEEDED ######
if [[ $writemass -eq 1 ]];then
echo '
2
Particle Data: AtomName(i) AtomicNumber(i) AtomicMass(i) GaussianWidth(i) 
O 8 16.0 15.0
H 1 1.  6.0
'>>FMSINPOUT/Control.dat
fi
############ END-OF-USER-INPUT ###########################

export TMPDIR=/scratch/$USER/MOLPROFMS_$$
$MOLPROEXE -s --no-xml-output -W $TMPDIR  >& $inputdeck.out < $inputdeck
rm -r /scratch/$USER/MOLPROFMS_$$

head -2 FMSINPOUT/TrajDump.1 | tail -1 > FMSTRAJS/Traj.$i
cp FMSINPOUT/FMS.out FMSTRAJS/FMS.out.$i 

#---The following not needed as long as MOLPRO executes without error
#---Must be in sync with make_restart.f90
#awk 'BEGIN{printf "%f ",0.0}{if($1=="Position:"){printf "%f ",$2;getline;printf "%f ",$1;getline;printf "%f ",$1}}'\
#    FMSTRAJS/FMS.out.$i > FMSTRAJS/Traj.$i
#awk '{if($1=="Momentum:"){printf "%f ",$2;getline;printf "%f ",$1;getline;printf "%f ",$1}}'\
#    FMSTRAJS/FMS.out.$i >> FMSTRAJS/Traj.$i
#awk 'BEGIN{printf " \n"}{}' FMSTRAJS/FMS.out.$i >> FMSTRAJS/Traj.$i

 echo -n "$i "
 let i++
done
echo "" 

