#!/bin/bash
cd MM

timestep=$1
ibead=$2
geom=../geom_mm.dat.$ibead

natom=`cat $geom | wc -l`
let natom1=natom+1

let lines=timestep*natom1+natom1
if [[ $natom -eq 3 ]];then
   head -$lines engrad_mmhigh.dat.all.ref | tail -$natom1 > ../engrad_mm.dat.$ibead
else 
   head -$lines engrad_mmlow.dat.all.ref | tail -$natom1 > ../engrad_mm.dat.$ibead
fi

