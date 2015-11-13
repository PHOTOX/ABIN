#!/bin/bash
ABINEXE="$PWD/$1 -x mini.dat"

cd TESTS 

function dif_files {
local status=0
local cont
local cont_no
cont_no=1
for file in $* 
do
   if [[ -e $file.ref ]];then
      diff $file $file.ref > $file.diff
      if [[ $? -ge 1 ]];then
         diff -y -W 200  $file $file.ref > $file.diff
         ../numdiff.py $file.diff
      fi
      if [[ $? -ge 1 && $cont_no -eq 1 ]];then
         status=1
         echo "File $file differ. Continue? [y/n]"
         while true 
         do
            read cont
            if [[ $cont = "n" || $cont = "no" ]];then
               echo "Exiting..."
               exit 1
            elif [[ $cont = "y" || $cont = "yes" ]];then
               echo "Continuing..."
               cont_no=0
               break
            else 
               echo "Please enter 'y' or 'n'"
            fi
         done
      fi
   fi
done
return $status
}

function makeref {
echo "Making new reference files."
for file in $* 
do
   if [[ -e $file.ref ]];then
      mv $file $file.ref
   fi
done
}

function clean {
rm -rf $* output
rm -f *.diff
if [[ -e "restart.xyz.0.ref" ]];then
   cp restart.xyz.0.ref restart.xyz
fi
}

err=0

files=( bkl.dat phase.dat coef.dat phaserest.dat phaserest.?? nacmrest.dat nacmrest.dat.?? minimize.dat geom.mini.xyz temper.dat r.dat vel.dat cv.dat cv_dcv.dat  dist.dat angles.dat dihedrals.dat geom.dat.??? geom_mm.dat.??? DYN/OUT* MM/OUT* state.dat stateall.dat ERROR debug.nacm dotprod.dat pop.dat prob.dat PES.dat energies.dat est_energy.dat movie.xyz movie_mini.xyz restart.xyz.old restart.xyz restart.xyz.?? restart.xyz.? )

#EULER should check wf_thresh conditions
# TODO: Make test_readme.txt, with specifications of every test
# TODO: by default, use mmwater as a potential instead of NAB
# Make tests for NAB, MPI and CP2K
if [[ $2 == "sh" ]];then
   folders=( SH_EULER SH_RK4 SH_BUTCHER SH_RK4_PHASE )
elif  [[ $2 = "all" || $2 = "clean" ]];then
   folders=( CMD GLE SH_EULER SH_RK4 SH_BUTCHER SH_RK4_PHASE SH_TDC PIGLE PIMD ABINITIO SHAKE HARMON MINI QMMM )
   if [[ $3 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=NAB
   fi
   if [[ $4 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=MPI
   fi
   if [[ $5 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=CP2K
   fi
   if [[ $6 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=FFTW
   fi
else
   folders=$2
fi
#folders=( SHAKE )

echo "Running tests in directories:"
echo ${folders[@]}

for dir in ${folders[@]}
do
   if [[ ! -e $dir ]];then
      echo "Directory $dir not found. Skipping...."
      continue
   else
      echo "Entering directory $dir"
   fi
   cd $dir 

   clean ${files[@]}

   if [[ $2 = "clean" ]];then
      echo "Cleaning files in directory $dir "
      cd ../
      continue
   fi

   $ABINEXE > output #-i input.in
   #for testing restart
   if [[ -e input.in2 ]];then
      $ABINEXE -i input.in2 >> output
   fi

   if [[ $7 = "makeref" ]];then

      makeref ${files[@]}

   else

      dif_files ${files[@]}

   fi

   if [[ $? -ne "0" ]];then
      err=1
   else
      echo "PASSED"
      echo "======================="
   fi
   cd ../
done

cd ../
echo " "

if [[ $err -ne "0" ]];then
   echo "Some tests DID NOT PASS."
else
   echo "All tests PASSED."
fi

exit $err


