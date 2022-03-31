#/bin/bash

rm -f EXIT ERROR plumed_*.dat abin.out* bck.0.plumed.out
if [[ "$1" = "clean" ]];then
   exit 0
fi

ABINEXE=$1
ABININ=input.in
ABINOUT=abin.out

# This is needed for compatibility with new Plumed 2.8
# See: https://github.com/plumed/plumed2/pull/757
export PLUMED_DP2CUTOFF_NOSTRETCH=1

$ABINEXE -i $ABININ -x mini.xyz > $ABINOUT
$ABINEXE -i ${ABININ}2  -x mini.xyz > ${ABINOUT}2

# Compatibility for Plumed 2.8
sed -i 's/stretched-gaussian/gaussian/' plumed_hills.dat
