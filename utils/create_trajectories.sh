#!/bin/bash

# ---------------------------------------------------------------------------------
# create_trajectories.sh - Generate and execute a set ABIN simulations
#
# Initial geometries (and optionally velocities) are taken sequentially from a XYZ trajectory file.

# The trajectories are executed and stored in $folder.

# Files needed in this folder:
#	$inputdir : template directory with ABIN input files (mainly input.in and r.abin)
# 	abin-randomint PRNG program for generating random seeds, should be in your $PATH.
#---------------------------------------------------------------------------------

# Exit if undefined variable is used
set -u

#### SETUP ####
# Path to a XYZ file with initial geometries
movie=coords.xyz
# Path to XYZ file with initial velocities (optional)
# veloc=vels.xyz

# Starting index for the initial geometries
isample=1
# End index
nsample=100

# Folder name where the trajectories will be created
folder=MY_MOLECULE_TRAJS

# Directory with the input files for ABIN (input.in et al)
inputdir=TEMPLATE-$folder

# File with ABIN input parameters, we need this path
# so we can inject random number seed into it.
abin_input=$inputdir/input.in

# Random seed to generate random seeds for individual trajectories
# Set to a negative number for a time-based random seed.
irandom0=156863189

# Specify path to launch script that is normally submitted to the queuing system.
# Comment out this line if you're running locally.
launch_script=$inputdir/r.abin
# If you don't provide a launch script,
# we need a (preferably absolute) path to abin executable
# abin_exe=/path/to/abin

# Comment out this line if you don't want to run calculations yet
# submit_command="qsub -cwd -V -q nq -cwd "
# If you don't use queing system (like SLURM), use the following line
# submit_command=bash

# Number of batch jobs to submit to queue
# set only if you have more trajectories than jobs
# jobs=20
########## END OF SETUP ##########


function files_exist {
  for f in "$@";do
    if [[ -n ${f-} && ! -f $f ]]; then
      echo "ERROR: File '$f' does not exist!"
      exit 1
    fi
  done
}

function folders_exist {
  for d in "$@";do
    if [[ -n ${f-} && ! -d $d ]]; then
      echo "ERROR: Directory '$d' does not exist!"
      exit 1
    fi
  done
}

folders_exist "$inputdir"
files_exist "$movie" "${veloc-}" "$abin_input" "${launch_script-}"

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

# Generate random number generator seeds for individual trajectories
echo "Generating $nsample random integers for random seeds"
echo "abin-randomint --seed $irandom0 --num $nsample > iran.dat"
if ! abin-randomint --seed $irandom0 --num $nsample > iran.dat
then
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
   head -$offset $movie | tail -$natom2 > $folder/TRAJ.$i/initial.xyz
   if [[ -n "${veloc-}" ]];then
      head -$offset "$veloc" | tail -$natom2 > $folder/TRAJ.$i/veloc.in
   fi

   ## Now prepare input.in and r.abin
   irandom=$(head -$i iran.dat |tail -1)

   # TODO: Validate this step
   sed -r "s/irandom *= *[0-9]+/irandom=$irandom/" $abin_input > $folder/TRAJ.$i/input.in 

   cat > $folder/TRAJ.$i/r.$folder.$i << EOF
#!/bin/bash
INPUTPARAM=input.in
INPUTGEOM=initial.xyz
EOF

   if [[ -n ${veloc-} ]];then
      echo "INPUTVELOC=veloc.in" >> $folder/TRAJ.$i/r.$folder.$i
   fi

   if [[ -n ${launch_script-} ]];then
      grep -v -e '/bin/bash' -e "INPUTPARAM=" -e "INPUTGEOM=" -e "INPUTVELOC=" $launch_script >> $folder/TRAJ.$i/r.$folder.$i
   else
      if [[ -n ${veloc-} ]]; then
         echo "$abin_exe -i input.in -x initial.xyz -v veloc.in > abin.out 2>&1" > $folder/TRAJ.$i/r.$folder.$i
      else
         echo "$abin_exe -i input.in -x initial.xyz > abin.out 2>&1" > $folder/TRAJ.$i/r.$folder.$i
      fi
   fi

   echo "(cd TRAJ.$i && bash r.$folder.$i)" >> $folder/$folder.$isample.$j.sh

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
    if [[ $submit_command = "bash" ]];then
        echo "Launching $j calculations locally"
        submit_command="nohup $submit_command"
    else
        echo "Submitting $j calculations with: $submit_command"
    fi
    while [[ $k -le $j ]]
    do
        if [[ -f $folder.$isample.$k.sh ]];then
            $submit_command $folder.$isample.$k.sh &
        fi
        (( k++ ))
    done
    # Wait for submit commands to finish (they should be fast!)
    if [[ $submit_command != "bash" ]]; then
        wait
    fi
fi
