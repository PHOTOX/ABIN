#!/bin/bash
set -euo pipefail

# Parameters are passed from Makefile.
if [[ $# -ne 8 ]]; then
  echo "ERROR: Incorrect number of parameters passed to $0"
  echo "Invoked as:"
  echo "$0 $@"
  exit 1
fi

ABINEXE="$PWD/bin/$1 -x mini.xyz"
ABINOUT="abin.out"
TESTS="$2"
MPI="$3"
FFTW="$4"
PLUMED="$5"
CP2K="$6"
TCPB="$7"
ACTION="$8"

# NOTE: For MPI tests, we rely on the fact that
# MPI_PATH is defined in make.vars and exported in Makefile.
# It's not needed if you use system-wide MPI installation
# or if mpirun is in your PATH.

if [[ $ACTION = "makeref" && $TESTS = "all" ]];then
  echo "ERROR: You should not call makeref on all tests at once."
  echo "Specify a concrete test which you want to modify, e.g."
  echo "make makeref TEST=CMD"
  exit 1
fi

if [[ $TCPB = "TRUE" ]];then
  export LD_LIBRARY_PATH=${TCPB_LIB}:${LD_LIBRARY_PATH-}
fi

cd $(dirname $0)
TESTDIR=$PWD

function diff_files {
  local return_status=0
  local error_code
  local ref_file
  local test_file
  # Compare test results with all existing reference files
  local reference_files=$(ls *.ref)
  if [[ -z $reference_files ]];then
    echo "ERROR: No reference files were found"
    return 1
  fi

  for ref_file in $reference_files
  do
    test_file=$(basename $ref_file .ref)
    if [[ ! -f $test_file ]];then
      # The output file does not exist.
      # Something went seriously wrong, ABIN probably crashed prematurely.
      # No need for further checks, exit NOW.
      echo "ERROR: Could not find output file \"$test_file\""
      return 1
    fi

    error_code=0
    diff -q $test_file $ref_file > /dev/null || error_code=$?
    if [[ $error_code -ne 0 ]];then
       # The reference file is different, but maybe it's just numerical noise?
       error_code=0
       diff -y -W 500 $test_file $ref_file | grep -e '|' -e '<' -e '>' > $test_file.diff

       ../numdiff.py $test_file.diff || error_code=$?

       if [[ $error_code -ne 0 ]];then
          # The changes were bigger that the thresholds specified in numdiff.py
          return_status=$error_code
       fi
    fi
  done
  return $return_status
}

# Update already existing reference files.
# Called by `make makeref TEST=TEST_FOLDER`
# If you're creating a completely new test,
# you need to create the reference files manually.
function makeref {
  local ref_file
  local test_file
  local reference_files=$(ls *.ref)
  echo "Making new reference files."
  if [[ -z $reference_files ]];then
    echo "ERROR: No reference files were found."
    exit 1
  fi
  for ref_file in $reference_files
  do
    test_file=$(basename $ref_file .ref)
    if [[ ! -f $test_file ]];then
      echo "ERROR: Could not find output file \"$test_file\""
      exit 1
    fi
    echo "mv $test_file $ref_file"
    mv $test_file $ref_file
  done
}

function clean {
  if [[ -f "test.sh" ]];then
    ./test.sh clean
  fi
  rm -rf $*
  rm -f *.diff
}

# List of all possible ABIN output files.
# Used by `make testclean` to cleanup test directories.
output_files=( *.dat *.out ERROR movie.xyz forces.xyz velocities.xyz
geom.dat.??? geom_mm.dat.??? geom.mini.xyz nacmrest.dat.?? hopgeom.*.xyz bck.*
restart_sh.bin restart_sh.bin.old restart_sh.bin.?? restart.xyz.old restart.xyz.? restart.xyz.?? restart.xyz )

# Run all tests
if [[ $TESTS = "all" ]];then
   folders=(INIT CMD NHC-GLOBAL SHAKE \
            SH_EULER SH_RK4 SH_BUTCHER SH_RK4_PHASE \
            SH_IGNORE SH_NACM_FAIL SH_S0S1 SH_ENERGY_DIFF SH_ENERGY_DRIFT \
            SH_BUTCHER_PHASE SH_SIMPLE_RESCALE SH_FRUSTRATED \
            LZ_SS LZ_ST LZ_ENE \
            PIMD ABINITIO ABINITIO-FAIL MTS \
            LANGEVIN QT QT2 PIGLE PIGLE2 GLE-CANONICAL \
            HARMON MORSE DOUBLEWELL SPLINE MM MINI QMMM \
            ANALYZE_EXT CMDLINE WATER_FAIL ERMD)

   if [[ $MPI = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=REMD
      let index++
      folders[index]=TERAPI
      let index++
      folders[index]=TERAPI-PIMD
      let index++
      folders[index]=TERAPI-PIMD-PARALLEL
      let index++
      folders[index]=TERAPI-REMD
      let index++
      folders[index]=TERAPI-FAILS
      let index++
      folders[index]=TERAPI-SH-S0
      let index++
      folders[index]=TERAPI-LZ
   else
      let index=${#folders[@]}+1
      folders[index]=WITHOUT_MPI
   fi

   if [[ $CP2K = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=CP2K
      if [[ $MPI = "TRUE" ]];then
        let index++
        folders[index]=CP2K_MPI
      fi
   else
      let index=${#folders[@]}+1
      folders[index]=WITHOUT_CP2K
   fi

   if [[ $FFTW = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=PIGLET
      let index++
      folders[index]=PIGLET2
      let index++
      folders[index]=PILE
   else
      let index=${#folders[@]}+1
      folders[index]=WITHOUT_FFTW
   fi

   if [[ $PLUMED = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=PLUMED
   else
      let index=${#folders[@]}+1
      folders[index]=WITHOUT_PLUMED
   fi

   if [[ $TCPB = "TRUE" ]];then
      let index=${#folders[@]}+1
      folders[index]=TCPB_FAIL
   else
      let index=${#folders[@]}+1
      folders[index]=WITHOUT_TCPB
   fi

else

   # Only one test selected, e.g. by running
   # make test TEST=CMD
   folders=${TESTS}

fi

echo "Running tests in directories:"
echo ${folders[@]}

global_error=0

for dir in ${folders[@]}
do
   if [[ ! -d $dir ]];then
      echo "Directory $dir not found. Exiting prematurely."
      exit 1
   fi
   echo "Entering directory $dir"
   cd $dir

   # Always clean the test directory before runnning the test.
   clean ${output_files[@]}

   # If we just want to clean the directories,
   # we skip the the actual test here
   if [[ $ACTION = "clean" ]];then
      echo "Cleaning files in directory $dir"
      cd $TESTDIR
      continue
   fi

   # For special cases such as REMD, we need a more complicated test setup.
   # If a file 'test.sh' is present in the test directory we will use it.
   if [[ -f "test.sh" ]];then

      # Redirection to dev/null apparently needed for CP2K tests.
      # Otherwise, STDIN is screwed up. I have no idea why.
      # http://stackoverflow.com/questions/1304600/read-error-0-resource-temporarily-unavailable
      # TODO: Figure out a different solution
      #./test.sh $ABINEXE 2> /dev/null
      ./test.sh $ABINEXE || true

   else
      if [[ -f "velocities.in" ]];then
         $ABINEXE -v "velocities.in" > $ABINOUT 2>&1 || true
      else
         $ABINEXE > $ABINOUT 2>&1 || true
      fi

      #for testing restart
      if [[ -e input.in2 ]];then
         $ABINEXE -i input.in2 >> $ABINOUT 2>&1 || true
      fi
   fi

   if [[ $ACTION = "makeref" ]];then

      makeref

   else

      # Since we're running in the -e mode,
      # we need to "hide" this possibly failing command
      # https://stackoverflow.com/a/11231970/3682277
      current_error=0
      diff_files || current_error=$?
      if [[ $current_error -ne 0 ]];then
        global_error=1
        echo -e "$dir \033[0;31mFAILED\033[0m"
      else
        echo -e "\033[0;32mPASSED\033[0m"
      fi
   fi

   echo "======================="

   cd $TESTDIR
done

echo " "

if [[ ${global_error} -ne 0 ]];then
   echo -e "Some tests \033[0;31mFAILED\033[0m."
else
   echo -e "\033[0;32mAll tests PASSED.\033[0m"
fi

exit $global_error
