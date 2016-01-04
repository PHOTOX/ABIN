#!/bin/bash

# A simle script for analysis of SH trajectories. 
# It Assumes that the trajectories are in folders $folder/TRAJ.$i
# It works even for trajectories that are still running by using fetchabin.sh.
# You should have the following scripts in your path: fetchabin.sh, checkenergy.awk, prum.awk
# It produces 3 files in folder $folder :
# energycons.dat - contains information about energy conservation in each trajectory
#                - 2nd column contains info about maximum energy difference in two consecutive timesteps
#                - 3rd column contains maximum energy drift of the trajectory
#                - 4th column contains final energy drift of the trajectory
# populations.dat - contains electronic populations based on actual state occupation from all trajectories.
# elpop.dat - contains electronic populations based on WF coefficients from file pop.dat
#           - It should be similar to populations.dat. If that's not the case, the SH algorithm does not work properly.

#######SETUP#############
folder=CAS22SA3.2
isample=1
nsample=100
SKIPFOLDERS=( 6 74 )     # these trajectories will be excluded from analysis
# The following variables are determined automatically from TRAJ.1
nstate=$(awk -F"[,=]" '{if($1=="nstate")print $2}' $folder/TRAJ.1/input.in)
dt=$(awk -F"[,=]" '{if($1=="dt")print $2}' $folder/TRAJ.1/input.in)
nstep=$(awk -F"[,=]" '{if($1=="nstep")print $2}' $folder/TRAJ.1/input.in)
##########END OF SETUP##########

shopt -s expand_aliases
alias grp="grep -s -v -e \"#\" " # used for filtering out the "#" characters

if [[ ! -d $folder ]];then
   echo "Folder $folder does not exist. Exiting..."
   exit 1
fi

if [[ -z $nstate ]] || [[ -z $dt ]] || [[ -z $nstep ]];then
   echo "Variable nstate, nstep or dt is not defined. Exiting."
   exit 1
fi

cd $folder
rm -f energycons.dat populations.dat elpop.dat pom.dat pom2.dat
cd ..

i=$isample

while [[ $i -le $nsample ]];do

   # skip trajectories specified by the user
   for j in ${SKIPFOLDERS[@]};do
      if [[ $j -eq $i ]];then
         echo "Omittind directory R$i"
         let i++
         continue
      fi
   done

echo -n "$i "
DEST=$folder/TRAJ.$i/

# Try to update with fetchabin.sh if the job is still running
if [[ -e $DEST ]];then
   cd $DEST
   grep -q -s "Job finished" output
   if [[ -e job.log && $? -ne 0 && ! -e ERROR ]];then
      fetchabin.sh  > /dev/null
   fi
   cd ../../
else
   echo "Folder $DEST does not exist! Ignoring..."
fi

# Check for the existence of pop.dat
grp -q $DEST/pop.dat 
if [[ $? -ne "0" ]];then
        echo "Problem reading pop.dat for TRAJ.$i. Skipping..."
	let i++
        continue
fi

###CHECKING ENERGY CONSERVATION
grp $DEST/energies.dat | checkenergy.awk -v itrj=$i >> $folder/pom2.dat

# Get populations from pop.dat
grp -e "deltaE" $DEST/pop.dat | \
   awk -v nst=$nstate '{printf"%i %i ",NR,$2;for(i=3;i<=nst+2;i++){printf"%f ",$i};printf"\n"}' >> $folder/pom.dat


let i++
done
echo  " "


# Analyzing energy conservation, getting averages
echo "# Traj.no. Max_energy_jump Max_energy_drift FinalDrift" > $folder/energycons.dat
cat $folder/pom2.dat >> $folder/energycons.dat
prum.awk -v column=2 $folder/pom2.dat >> $folder/energycons.dat
prum.awk -v column=3 $folder/pom2.dat >> $folder/energycons.dat


# ANALYZING POPULATIONS FROM ALL TRAJECTORIES

awk -v dt="$dt" -v nstate="$nstate" -v ngeom="$nstep" -v file="$folder/elpop.dat" 'BEGIN{
for(i=1;i<nstate;i++) {for(j=1;j<=ngeom;j++) pop[i,j]=0.0}}
{
   #Populations based on state occupation
   pop[$2,$1]++
   #populations based on el. coefficients
   for (i=3;i<=nstate+2;i++)
      elpop[i-2,$1]=elpop[i-2,$1]+$i
}

END{
print "# Time [fs] Populations based on state occupations"
print "# Time [fs] Populations based on electronic coefficients" >> file
print "# #Number of living trajectories at each timestep"
for(i=1;i<=ngeom;i++) {
   anorm=0.0
   for (j=1;j<=nstate;j++) anorm+=pop[j,i]
   if (anorm<0.5) break
   print "#",anorm
   printf"%f ", i*dt*0.02419
   printf"%f ", i*dt*0.02419 >> file
   for (j=1;j<=nstate;j++) {
      printf" %f",pop[j,i]/anorm
   }
   printf"\n"

   #printinf electronic populations to a different file
   #we are assuming anorm is the same, i.e. sum of populations is near 1 all the time
   for (j=1;j<=nstate;j++) {
      printf" %f",elpop[j,i]/anorm >> file
   }
   printf"\n" >> file
}
}' $folder/pom.dat > $folder/populations.dat

