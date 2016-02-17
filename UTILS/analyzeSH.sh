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
folder=TERACAS-22-SA3-MPI
isample=1
nsample=100
nstate=2
dt=5.0
nstep=1700
SKIPFOLDERS=(   )     # these trajectories will be excluded from analysis
##########END OF SETUP##########
outfolder=$folder/ANALYSIS

currfolder=$PWD

echo "Number of steps: $nstep "
echo "Time step: $dt" 
echo "Number of states: $nstate"

shopt -s expand_aliases
alias grp="grep -s -v -e \"#\" " # used for filtering out the "#" characters

if [[ ! -d $folder ]];then
   echo "Folder $folder does not exist. Exiting..."
   exit 1
fi

if [[ ! -d $outfolder ]];then
   mkdir $outfolder
fi

if [[ -z $nstate ]] || [[ -z $dt ]] || [[ -z $nstep ]];then
   echo "Variable nstate, nstep or dt is not defined. Exiting."
   exit 1
fi

cd $outfolder
rm -f alive.dat simulation_times.dat energycons.dat populations.dat elpop.dat pom.dat pom2.dat
cd $currfolder

i=$isample

while [[ $i -le $nsample ]];do

   # skip trajectories specified by the user
   for j in ${SKIPFOLDERS[@]};do
      if [[ $j -eq $i ]];then
         echo "Omitting directory R$i"
         let i++
         continue
      fi
   done

   echo -n "$i "
   DEST=$folder/TRAJ.$i/

   # Try to update with fetchabin.sh if the job is still running
   if [[ -d $DEST ]];then
      cd $DEST
      grep -q -s "Job finished" output
      if [[ -e job.log && $? -ne 0 && ! -e ERROR ]];then
         fetchabin.sh  > /dev/null
      fi
      cd $currfolder
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
   grp $DEST/energies.dat | checkenergy.awk -v itrj=$i >> $outfolder/pom2.dat
   
   # Get populations from pop.dat
   grp -e "deltaE" $DEST/pop.dat | \
   awk -v nst=$nstate '{printf"%i %i ",NR,$2;for(i=3;i<=nst+2;i++){printf"%f ",$i};printf"\n"}' >> $outfolder/pom.dat

#  Checking total simulation time in fs
   grp $DEST/energies.dat | tail -n 1 | awk  -v "i=$i" '{print i, $1}' >> $outfolder/simulation_times.dat

   let i++
done
echo  " "


# Analyzing energy conservation, getting averages
echo "# Traj.no. Max_energy_jump Max_energy_drift FinalDrift" > $outfolder/energycons.dat
cat $outfolder/pom2.dat >> $outfolder/energycons.dat
prum.awk -v column=2 $outfolder/pom2.dat >> $outfolder/energycons.dat
prum.awk -v column=3 $outfolder/pom2.dat >> $outfolder/energycons.dat


# ANALYZING POPULATIONS FROM ALL TRAJECTORIES
# Estimate standard error for populations based on binomial distributions
# See https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
# and http://stats.stackexchange.com/questions/111355/confidence-interval-and-sample-size-multinomial-probabilities
# http://sites.stat.psu.edu/~sesa/stat504/Lecture/lec3_4up.pdf

echo $nstep 
awk -v dt="$dt" -v nstate="$nstate" -v ngeom="$nstep" -v file="$outfolder/elpop.dat" -v file2="$outfolder/alive.dat" 'BEGIN{
for(i=1;i<nstate;i++) {
   for(j=1;j<=ngeom;j++) {
      pop[i,j]=0.0
      elpop[i,j]=0.0
      elpop2[i,j]=0.0
   }
}
}

{
   #Populations based on state occupation
   pop[$2,$1]++
   #populations based on el. coefficients
   for (i=3;i<=nstate+2;i++) {
      elpop[i-2,$1]+=$i
      elpop2[i-2,$1]+=$i*$i
   }
}

END{
print "# Time [fs] No._Trajs  Populations based on state occupations"
print "# Time [fs] Populations based on electronic coefficients" >> file
print "# Time [fs]  Number of living trajectories" >> file2

# NORMALIZATION
for(i=1;i<=ngeom;i++) {
   anorm[i]=0.0

   for (j=1;j<=nstate;j++)  anorm[i]+=pop[j,i]

   if (anorm[i]<0.5) break

   for (j=1;j<=nstate;j++) {
      pop[j,i]/=anorm[i]
      elpop[j,i]/=anorm[i]
   }
   print i*dt*0.02419, anorm[i] >> file2
}

# PRINTING
for (j=1;j<=nstate;j++) {
   for(i=1;i<=ngeom;i++) {

      print i*dt*0.02419,pop[j,i],sqrt(pop[j,i]*(1-pop[j,i])/anorm[i])*1.96

      # printing electronic populations to a different file
      # we are assuming anorm is the same, i.e. sum of populations is near 1 all the time
      print i*dt*0.02419,elpop[j,i],sqrt( elpop2[j,i]/anorm[i]-(elpop[j,i])^2) >> file
   }
   print ""
   print ""
   print "" >> file
   print "" >> file
}
}' $outfolder/pom.dat > $outfolder/populations.dat

cd $outfolder
# Make a plot, using also Standard errors
cat > gnuplot.in << EOF
set term png enhanced size 800,600
set output "populations.png"
unset key
set multiplot 
set size 1.0, 1.0 
set origin 0.0, 0.0
set tics nomirror
totaltime=$nstep*$dt*0.02419

#TODO: make color styles
# same color for lines and points with error bars

set yrange [0:1]
set xrange [0:totaltime]
set ylabel "Populations [ - ]" 
set xlabel "Time [ fs ]" 
EOF
for((i=0;i<nstate;i++));do
cat >> gnuplot.in << EOF
plot "populations.dat" index $i every 100 using 1:2:3 with yerrorbars,\
     "populations.dat" index $i using 1:2 with lines
EOF
done

gnuplot  gnuplot.in 


cd $currfolder
