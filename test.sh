#!/bin/bash
ABINEXE="$PWD/$1 -x mini.dat"

function dif_files {
local status=0
local cont
local cont_no
cont_no=1
for file in $* 
do
   if [[ -e $file.ref ]];then
      diff  $file $file.ref > $file.diff
      if [[ $? -ge 1 && $cont_no -eq 1 ]];then
         status=1
         echo "File $file differ. Continue? [y/n]"
         while true 
         do
            read cont
            if [[ $cont = "n" ]];then
               echo "Exiting..."
               exit 1
            elif [[ $cont = "y" ]];then
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
rm -f $* output
rm -f *.diff
if [[ -e "restart.xyz.0.ref" ]];then
   cp restart.xyz.0.ref restart.xyz
fi
}


cd TESTS
err=0

files=( bkl.dat phase.dat coef.dat phaserest.dat phaserest.?? nacmrest.dat nacmrest.dat.?? minimize.dat geom.mini.xyz temper.dat r.dat vel.dat cv.dat cv_dcv.dat  dist.dat angles.dat dihedrals.dat geom.dat.??? geom_mm.dat.??? ORCA/input* MM/input* state.dat stateall.dat ERROR debug.nacm dotprod.dat pop.dat prob.dat PES.dat energies.dat est_energy.dat movie.xyz movie_mini.xyz restart.xyz.old restart.xyz restart.xyz.?? restart.xyz.? )

#EULER should check wf_thresh conditions
if [[ $2 == "sh" ]];then
folders=( SH_EULER SH_RK4 SH_BUTCHER SH_RK4_PHASE )
else
folders=( CMD SH_EULER SH_RK4 SH_BUTCHER SH_RK4_PHASE SH_TDC PIGLE GLE PIMD ABINITIO SHAKE HARMON MINI QMMM )
#folders=( SHAKE )
fi

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
      cd ../
      continue
   fi

   $ABINEXE > output #-i input.in
   #for testing restart
   if [[ -e input.in2 ]];then
      $ABINEXE -i input.in2 >> output
   fi

   if [[ $2 = "makeref" ]];then

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


