# !/bin/bash

# This is a sample script which takes
# geometries consecutively from xyz trajectory and launches G09 jobs.

# You should use PickGeoms.sh before running this script.

# It is advisable to use the timestep from the movies as a comment.

#########SETUP########
name=wat             # name of the job
nproc=2              # how many processors should we use?
nsample=20            # last geometry
ifirst=1             # first geometry
movie=geoms.xyz      # file with xyz geometries
spin=1               # molecular spin
charge=0             # molecular charge
nrun=8              # determines number of jobs to submit
#submit="qsub -q sq-8-16"
######################


if [[ ! -e $movie ]];then
   echo "ERROR: File $movie does not exist."
   exit 1
fi

if [[ $nrun -gt $nsample ]];then
   echo "WARNING: Number of jobs is bigger than number of samples."
   nrun=$nsample
fi

natom=$(head -1 $movie )            # number of atoms
let natom2=natom+2
let natom1=natom+1

lines=`cat $movie | wc -l` 
geoms=`expr $lines / $natom2`

# determine number of G09 calculations per job
let njobs=nsample/nrun
#determine the remainder and distribute it evenly between jobs
let remainder=nsample-njobs*nrun

if [[ $nsample -gt $geoms ]];then
   echo "ERROR: Number of geometries is smaller than number of samples."
   echo "Change parameter \"nsample\"."
   exit 1
fi

rm -f r.$name.*

j=1
offset=0
i=$ifirst

while [[ $i -le $nsample ]]
do
   let offset=offset+natom2
# G09 STUFF- MODIFY TO YOUR NEEDS
cat > $name.$i.com <<EOF
\$rungauss
%Mem=520Mb
%NProcShared=$nproc
#BHandHLYP/6-31++g** sp

EOF

#--Put time step in the comment
   head -$offset $movie | tail -$natom1 | head -1 >> $name.$i.com

   echo " " >> $name.$i.com
   echo $charge $spin >> $name.$i.com

   head -$offset $movie | tail -n $natom >> $name.$i.com
   echo " " >>$name.$i.com


   echo "~hollas/bin/G09 $name.$i.com" >>r.$name.$j
#----END OF G09 STUFF


#--Distribute G09 calculations evenly between jobs for queue
   if [[ $remainder -le 0 ]];then
      let ncalc=njobs
   else
      let ncalc=njobs+1 
   fi
   if [[ `expr $i % $ncalc` -eq 0 ]] && [[ $j -lt $nrun ]]; then
      let j++
      let remainder--
   fi

   let i++

done

j=1

# SUBMIT JOBS
if [[ ! -z $submit ]];then
   while [[ $j -le $nrun ]]
   do
      $submit -cwd -pe shm $nproc r.$name.$j
      let j++
   done
fi

