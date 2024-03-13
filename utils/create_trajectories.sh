#!/bin/bash

# ---------------------------------------------------------------------------------
# create_trajectories.sh - Generate and execute a set ABIN simulations
#
# Initial geometries (and optionally velocities) are taken sequentially from XYZ movie files.

# The trajectories are executed and stored in $folder.

# Files needed in this folder:
#	$inputdir : template directory with ABIN input files (mainly input.in and r.abin)
# 	abin-randomint PRNG program for generating random seeds, should be in your $PATH.
#---------------------------------------------------------------------------------

# Exit if undefined variable is used
set -u

#######-----SETUP---#############
irandom0=156863189  # random seed, set negative for random seed based on time
movie=coords.xyz    # PATH TO a XYZ movie with initial geometries
veloc=vels.xyz      # PATH to XYZ initial velocities (optional)
isample=1	    # initial number of traj
nsample=100	    # number of trajectories
folder=MP2-NH4      # Name of the folder with trajectories
inputdir=TEMPLATE-$folder   # Directory with input files for ABIN
abin_input=$inputdir/input.in   # main input file for ABIN
launch_script=$inputdir/r.abin	# this is the file that is submitted by qsub

# Comment out this line if you don't want to run calculations yet
# submit_command="qsub -cwd -V -q nq -cwd "
# If you don't use queing system (like SLURM), use the following line
# submit_command=bash

# Number of batch jobs to submit, set only if you have more trajectories than jobs
# jobs=20
########## END OF SETUP ##########


function files_exist {
  for f in "$@";do
    if [[ ! -f $f ]];then
      echo "ERROR: File '$f' does not exist!"
      exit 1
    fi
  done
}

function folders_exist {
  for d in "$@";do
    if [[ ! -d $d ]];then
      echo "ERROR: Directory '$d' does not exist!"
      exit 1
    fi
  done
}

folders_exist "$inputdir"
files_exist "$movie" "$abin_input" "$launch_script"
if [[ -n "$veloc" ]];then
  files_exist "$veloc"
fi

natom=$(head -1 $movie)
if [[ $natom -lt 1 ]];then
  echo "ERROR: Invalid number of atoms on the first line of file $movie"
  exit 1
fi
echo "Number of atoms = $natom"

# TODO: Verify number of atoms and lines in the velocity file
(( natom2=natom+2 ))
lines=$(wc -l < $movie)
(( geoms=lines/natom2 ))
if [[ $nsample -gt $geoms ]];then
   echo "ERROR: Number of geometries ($geoms) is smaller than number of samples($nsample)."
   echo "Change parameter \"nsample\"."
   exit 1
fi

# determine number of ABIN simulations per job
(( nsimul=nsample-isample+1 ))
if [[ -z ${jobs-} ]]; then
   jobs=$nsimul
fi

if [[ $nsimul -le $jobs ]];then
   remainder=0
   injob=1
   jobs=$nsimul
else
   (( injob=nsimul/jobs ))  #number of simulations per job
   # determine the remainder and distribute it evenly between jobs
   (( remainder=nsimul-injob*jobs ))
fi

j=1
i=$isample
w=0 #current number of simulations in current j-th job

#--------------------generation of random numbers--------------------------------
echo "Generating $nsample random integers for random seeds"
echo "abin-randomint --seed $irandom0 --num $nsample > iran.dat"
abin-randomint --seed $irandom0 --num $nsample > iran.dat
if [[ $? -ne "0" ]];then
  echo "ERROR: Could not generate random numbers"
  exit 1
fi

mkdir -p $folder
cp iseed0 "$abin_input" $folder

(( offset=natom2*isample-natom2 ))

while [[ $i -le "$nsample" ]];do

   (( offset=offset+natom2 ))

   if [[ -d "$folder/TRAJ.$i" ]];then

         echo "Trajectory number $i already exists!"
         echo "If you want to overwrite it, first remove it:"
         echo "'rm -r $folder/TRAJ.$i'"
         exit 1

   else

      mkdir $folder/TRAJ.$i

   fi

   # Copy all the files from the template directory
   cp -r $inputdir/* $folder/TRAJ.$i

   # Prepare input geometry and velocities
   head -$offset $movie | tail -$natom2 > initial.xyz
   if [[ -n "${veloc-}" ]];then
      head -$offset $veloc | tail -$natom2 > $folder/TRAJ.$i/veloc.in
   fi

   ## Now prepare input.in and r.abin
   irandom=$(head -$i iran.dat |tail -1)

   # TODO: Validate this step
   sed -r "s/irandom *= *[0-9]+/irandom=$irandom/" $abin_input > $folder/TRAJ.$i/input.in 

   cat > $folder/TRAJ.$i/r.$folder.$i << EOF
#!/bin/bash
JOBNAME=ABIN.$folder.${i}_$$_\${JOB_ID}
INPUTPARAM=input.in
INPUTGEOM=initial.xyz
OUTPUT=abin.out
EOF

   if [[ -n ${veloc-} ]];then
      echo "INPUTVELOC=veloc.in" >> $folder/TRAJ.$i/r.$folder.$i
   fi

   grep -v -e '/bin/bash' -e "JOBNAME=" -e "INPUTPARAM=" -e "INPUTGEOM=" -e "INPUTVELOC=" $launch_script >> $folder/TRAJ.$i/r.$folder.$i

   chmod 755 $folder/TRAJ.$i/r.$folder.$i


   echo "cd TRAJ.$i || exit" >> $folder/$folder.$isample.$j.sh
   echo "./r.$folder.$i" >> $folder/$folder.$isample.$j.sh
   echo "cd $PWD/$folder || exit" >> $folder/$folder.$isample.$j.sh

   # Distribute calculations evenly between jobs for queue
   if [[ $remainder -le 0 ]];then
      ncalc=injob
   else
      (( ncalc=injob+1 ))
   fi
   (( w++ ))
   if [[ $w -eq $ncalc ]] && [[ $j -lt $jobs ]]; then
      w=0
      (( j++ ))
      (( remainder-- ))
   fi

   (( i++ ))

done

# Submit jobs
k=1
if [[ -n "${submit_command-}" ]];then
   cd $folder || exit 1
   while [[ $k -le $j ]]
   do
      if [[ -f $folder.$isample.$k.sh ]];then
         $submit_command $folder.$isample.$k.sh
      fi
      (( k++ ))
   done
fi
