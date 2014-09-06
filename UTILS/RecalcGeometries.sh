#!/bin/bash

# This is a sample script which takes geometries
# consecutively from xyz trajectory and launches G09 jobs.

# You should use PickGeoms.sh before running this script.

# EXAMPLE: if you have 10 geometries, and you want to skip first 3 of them
# and calculate the rest, set "first=4" and "nsample=0"
# in this case, "nsample=0" will do the same as "nsample=7"

# If you want to use some other program, you have to modify
# the part of the script beginning with "G09 STUFF".

#########SETUP########
name=test            # name of the job
first=3              # first geometry, will skip (first-1)geometries
nsample=5            # number of geometries, positive integer or 0 for all geometries from the one set as first
movie=geoms.xyz      # file with xyz geometries
g09="#BMK/aug-cc-pVDZ gfinput IOP(6/7=3) nosymm TD=(singlets,nstate=5)"
spin=1               # molecular spin
charge=0             # molecular charge
nproc=1              # how many processors should we use?
mem=500Mb            # memory in G09 job
jobs=4               # determines number of jobs to submit to queue
#submit="qsub -q sq-8-16"  # comment this line if you do not want to submit jobs
######################

g09path="/home/hollas/bin/G09"

if [[ ! -e $movie ]];then
   echo "ERROR: File $movie does not exist."
   exit 1
fi

natom=$(head -1 $movie )            # number of atoms
let natom2=natom+2
let natom1=natom+1

lines=`cat $movie | wc -l` 
geoms=`expr $lines / $natom2`

if [[ $nsample -eq 0 ]];then
   let nsample=geoms-first+1
fi

if [[ $jobs -gt $nsample ]];then
   echo "WARNING: Number of jobs is bigger than number of samples."
   jobs=$nsample
fi

last=`expr $first + $nsample - 1`

# determine number of G09 calculations per job
let injob=nsample/jobs
#determine the remainder and distribute it evenly between jobs
let remainder=nsample-injob*jobs

if [[ $nsample -gt `expr $geoms - $first + 1` ]];then
   echo "ERROR: Number of geometries ($geoms) decreased by unused geometries at the beginning is smaller than  number of samples."
   echo "Change parameter \"nsample\" or \"first\"."
   exit 1
fi

rm -f r.$name.*

j=1
let offset=(first-1)*natom2
i=$first

########################################################################
# YOU WILL NEED TO CHANGE THE FOLLOWING IF YOU WANT TO USE OTHER PROGRAM

# G09 STUFF - MODIFY TO YOUR NEEDS
cat > $name.template.g09 <<EOF
%Mem=$mem
%NProcShared=$nproc
$g09

EOF

while [[ $i -le $last ]]
do
   let offset=offset+natom2

   cp $name.template.g09 $name.$i.com

# It is advisable to use the timestep from the movies as a comment for future reference.
#-Put time step in the comment
   head -$offset $movie | tail -$natom1 | head -1 >> $name.$i.com

   echo " " >> $name.$i.com
   echo $charge $spin >> $name.$i.com

   head -$offset $movie | tail -n $natom >> $name.$i.com
   echo " " >>$name.$i.com


   echo "$g09path $name.$i.com" >>r.$name.$j
#----END OF G09 STUFF



########################################################################
##  THE REST IS GENERAL AND DOES NOT NEED TO BE MODIFIED


#--Distribute G09 calculations evenly between jobs for queue
   if [[ $remainder -le 0 ]];then
      let ncalc=injob
   else
      let ncalc=injob+1 
   fi
   if [[ `expr \( $i - $first + 1 \) % $ncalc` -eq 0 ]] && [[ $j -lt $jobs ]]; then
      let j++
      let remainder--
   fi

   let i++

done

j=1

# SUBMIT JOBS
if [[ ! -z $submit ]];then
   while [[ $j -le $jobs ]]
   do
      $submit -cwd -pe shm $nproc r.$name.$j
      let j++
   done
fi

