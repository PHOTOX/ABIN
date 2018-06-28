#!/bin/bash
ABINEXE="$PWD/$1 -x mini.dat"


MPI=$(awk -F"[# ,=]+" '{if($1=="MPI")print $2}' make.vars)
if [[ $MPI = "TRUE" ]];then
   export MPI_PATH=$(awk -F"[# ,=]+" '{if($1=="MPI_PATH") print $2}' make.vars)
   export LD_LIBRARY_PATH=$MPI_PATH/lib:$LD_LIBRARY_PATH
fi

cd TESTS 

function dif_files {
local status=0
local cont
local cont_no
local files
local f
cont_no=1
# Do comparison for all existing reference files
files=$(ls *.ref)
for f in $files   # $* 
do
   file=$(basename $f .ref)
   if [[ -e $file.ref ]];then  # this should now be always true
      diff $file $file.ref > $file.diff
      if [[ $? -ge 1 ]];then
         diff -y -W 320  $file $file.ref | egrep -e '|' -e '<' -e '>' > $file.diff
         ../numdiff.py $file.diff 2> /dev/null
      fi
      if [[ $? -ge 1 && $cont_no -eq 1 ]];then
         status=1
         echo "File $file differ from the reference. Continue anyway? [y/n]"
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
local files
local f
echo "Making new reference files."
files=$(ls *.ref)
for f in $files 
do
   file=$(basename $f .ref)
   if [[ -f $file.ref ]];then
      mv $file $file.ref
   else
      echo "Something horrible happened during makeref"
      exit 1
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

files=( WATER-RESTART.wfn* cp2k.out bkl.dat phase.dat wfcoef.dat phaserest.dat phaserest.?? nacmrest.dat nacmrest.dat.?? nacm_all.dat minimize.dat geom.mini.xyz temper.dat temper.dat radius.dat vel.dat cv.dat cv_dcv.dat  dist.dat angles.dat dihedrals.dat geom.dat.??? geom_mm.dat.??? DYN/OUT* MM/OUT* state.dat stateall.dat stateall_grad.dat ERROR debug.nacm dotprod.dat pop.dat prob.dat PES.dat energies.dat est_energy.dat movie.xyz movie.xyz movie_mini.xyz restart.xyz.old restart.xyz.? restart.xyz.?? restart.xyz )

# EULER should check wf_thresh conditions
# TODO: Make test_readme.txt, with specifications of every test, maybe include only in input.in
if [[ $2 == "sh" ]];then
   folders=( SH_EULER SH_RK4 SH_BUTCHER SH_RK4_PHASE )
elif  [[ $2 = "all" || $2 = "clean" ]];then
   folders=( CMD GLE SH_EULER SH_RK4 SH_BUTCHER SH_RK4_PHASE SH_TDC PIGLE PIMD ABINITIO SHAKE HARMON MINI QMMM )
   if [[ $3 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=NAB
      let index++
      folders[index]=NAB_HESS
   fi
   if [[ $4 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=REMD
#      let index++
      #      folders[index]=TERAPI # does not yet work
   fi
   if [[ $5 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=CP2K
      if [[ $4 = "TRUE" ]];then
         let index++
         folders[index]=CP2K_MPI
      fi
   fi
   if [[ $6 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=PIGLET
      let index++
      folders[index]=PILE
   fi
   if [[ $7 = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=PLUMED
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

   if [[ -f "test.sh" ]];then
      ./test.sh clean 
   else
      clean ${files[@]}
   fi

   if [[ $2 = "clean" ]];then
      echo "Cleaning files in directory $dir "
      cd ../
      continue
   fi

   # Redirection to dev/null apparently needed for CP2K tests.
   # Otherwise, STDIN is screwed up. I have no idea why.
   # http://stackoverflow.com/questions/1304600/read-error-0-resource-temporarily-unavailable
   if [[ -f "test.sh" ]];then

      ./test.sh $ABINEXE 2> /dev/null

   else
      if [[ -f "veloc.in.ref" ]];then
         $ABINEXE -v "veloc.in.ref" > output  
      else
         $ABINEXE > output  
      fi

      #for testing restart
      if [[ -e input.in2 ]];then
         $ABINEXE -i input.in2 >> output
      fi
   fi

   if [[ $8 = "makeref" ]];then

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


