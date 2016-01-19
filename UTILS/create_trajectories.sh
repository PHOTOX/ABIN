#!/bin/bash

#---------------------------------------------------------------------------------
#  Create_Trajectories                   Daniel Hollas, Ondrej Svoboda 2016

# This script generates and executes a set of dynamical trajectories using ABIN.
# Initial geometries are taken sequentially from a XYZ movie file.
# Initial velocities are optiononaly taken sequentially from a XYZ file.

# The script is designed both for surface hopping and adiabatic AIMD.

# The trajectories are executed and stored in $folder.

# Files needed in this folder:
#	$INPUTDIR : template directory with ABIN input files (mainly input.in and r.abin)
# 	traj-randomint PRNG program should be in your $PATH.
#---------------------------------------------------------------------------------

#######-----SETUP---#############
irandom0=-1          # random seed, set negative for random seed based on time
movie=geoms.xyz      # PATH TO a XYZ movie with initial geometries
veloc=velocs.dat     # PATH to XYZ initial velocities, leave blank if you do not have them
isample=1	     # initial number of traj
nsample=100	     # number of trajectories
INPUTDIR=TRAJS-TEMPLATE-QCHEM   # Directory with input files for ABIN
abin_input=$INPUTDIR/input.in   # main input file for ABIN
launch_script=$INPUTDIR/r.abin	# this is the file that is submitted by qsub


folder=MP2           # Name of the folder with trajectories
abinexe=$PHOTOX/bin/abin      # path to abin binary
submit="qsub -q nq -cwd  " # -q sq-* "    # comment this line if you don't want to submit to queue yet
rewrite=0            # if =1 -> rewrite trajectories that already exist
jobs=30              # number of batch jobs to submit. Trajectories will be distributed accordingly.

# Number of atoms is determined automatically from input.in
natom=$(awk -F"[! ,=]+" '{if($1=="natom")print $2}' $abin_input) #number of atoms
molname=$folder # Name of the job in the queue

## If you need to process the input via cut_sphere, adjust the following command
# cut="cut_sphere -u 4 -v 3" # cut command, velocity file is handled automatically if needed
# (Do not provide file names.)
##########END OF SETUP##########


function Folder_not_found {
   echo "Error: Folder $1 does not exists!"
   exit 1
}

function File_not_found {
   echo "Error: File $1 does not exists!"
   exit 1
}

function Error {
   echo "Error from command $1. Exiting!"
   exit 1
}

if [[ ! -d "$INPUTDIR" ]];then
   Folder_not_found "$INPUTDIR"
fi

if [[ ! -f "$movie" ]];then
   File_not_found "$movie"
fi

if [[ ! -z "$veloc" ]] && [[ ! -f "$veloc" ]];then
   File_not_found "$veloc"
fi

if [[ ! -e "$abinexe" ]];then
   File_not_found "$abinexe"
fi

if [[ ! -e $abin_input ]];then
   File_not_found "$abin_input"
fi

if [[ ! -e "$launch_script" ]];then
   File_not_found "$launch_script"
fi

if [[ -e "mini.dat" ]] || [[ -e "restart.xyz" ]];then
   echo "Error: Files mini.dat and/or restart.xyz were found here."
   echo "Please remove them."
   exit 1
fi

#   -------------------------------------------------------------------------------------

echo "Number of atoms = $natom"

# I don't think this is needed, we make sure not to overwrite any trajectory anyway
#if [[ -e $folder/$molname.$isample.*.sh ]];then
#   echo  "Error: File $folder/$molname.$isample.*.sh already exists!"
#   echo "Please, make sure that it is not currently running and delete it."
#   exit 1
#fi

# determine number of ABIN simulations per job
let nsimul=nsample-isample+1
if [[ $nsimul -le $jobs ]];then
   remainder=0
   injob=1
   jobs=$nsimul
else
   let injob=nsimul/jobs  #number of simulations per job
   # determine the remainder and distribute it evenly between jobs
   let remainder=nsimul-injob*jobs
fi


j=1
i=$isample

let natom2=natom+2

#--------------------generation of random numbers--------------------------------
echo "Generating $nsample random integers."
trajs-randomint $irandom0 $nsample > iran.dat
if [[ $? -ne "0" ]];then
   Error "trajs-randomint"
fi

#--------------------------------------------------------------------------------

mkdir -p $folder
cp iseed0 "$abin_input" $folder

let offset=natom2*isample-natom2

while [[ $i -le "$nsample" ]];do

   let offset=offset+natom2     

   if [[ -d "$folder/TRAJ.$i" ]];then
      if [[ "$rewrite" -eq "1" ]];then

         rm -r $folder/TRAJ.$i ; mkdir $folder/TRAJ.$i

      else

         echo "Trajectory number $i already exists!"
         echo "Exiting..."
         exit 1

      fi

   else

      mkdir $folder/TRAJ.$i

   fi

   cp -r $INPUTDIR/* $folder/TRAJ.$i


#--Now prepare mini.dat (and possibly veloc.in)

   head -$offset $movie | tail -$natom2 > geom
   if [[ ! -z "$veloc" ]];then
      head -$offset $veloc | tail -$natom2 > veloc.in
   fi
#--- Perform cutting, if needed.
   if [[  ! -z "$cut" ]];then
      if [[ ! -z "$veloc" ]];then
         $cut -vel veloc.in < geom
      else
         $cut
      fi
      if [[ $? -ne "0" ]];then
         Error "$cut"
      fi

      cp cut_qm.xyz geom
      cp veloc.out veloc.in
   fi

   mv geom $folder/TRAJ.$i/mini.dat

   if [[ ! -z "$veloc" ]];then
      mv veloc.in $folder/TRAJ.$i/
   fi


## Now prepare input.in and r.abin
   irandom=`head -$i iran.dat |tail -1`

   sed -r "s/irandom *= *[0-9]+/irandom=$irandom/" $abin_input > $folder/TRAJ.$i/input.in 

   cat > $folder/TRAJ.$i/r.$molname.$i << EOF
#!/bin/bash
ABINEXE=$abinexe
JOBNAME=ABIN.$molname.${i}_$$_\${JOB_ID}
INPUTPARAM=input.in
INPUTGEOM=mini.dat
OUTPUT=output
EOF

   if [[ ! -z $veloc ]];then
      echo "INPUTVELOC=veloc.in" >> $folder/TRAJ.$i/r.$molname.$i
   fi

   grep -v -e '/bin/bash' -e 'ABINEXE=' -e "JOBNAME=" -e "INPUTPARAM=" -e "INPUTGEOM=" -e "INPUTVELOC=" $launch_script >> $folder/TRAJ.$i/r.$molname.$i

   chmod 755 $folder/TRAJ.$i/r.$molname.$i


   echo "cd TRAJ.$i" >> $folder/$molname.$isample.$j.sh
   echo "./r.$molname.$i" >> $folder/$molname.$isample.$j.sh
   echo "cd $PWD/$folder" >> $folder/$molname.$isample.$j.sh

#--Distribute calculations evenly between jobs for queue
   if [[ $remainder -le 0 ]];then
      let ncalc=injob
   else
      let ncalc=injob+1 
   fi
   if [[ `expr \( $i - $first + 1 \) % $ncalc` -eq 0 ]] && [[ $j -lt $jobs ]]; then
      let j++
      let remainder--
   fi
#---------------------------------------------------------------------------

   let i++

done

# Submit jobs

if [[ ! -z "$submit" ]];then
   cd $folder
   while [[ $k -le $j ]]
   do
      if [[ -f $molname.$isample.$k.sh ]];then
         $submit -V -cwd $molname.$isample.$k.sh
      fi
      let k++
   done
fi

