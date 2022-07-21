#!/bin/bash

# This script must reside in folder 'MM'!

cd $(dirname $0)
source ../SetEnvironment.sh OPENMM

timestep=$1
ibead=$2
input=input$ibead
geom=../geom_mm.dat.$ibead
natom=$(wc -l < $geom)

# TBD: create PDB file, probably according to some template
templatefile=template.pdb

# Also, we need to take care of the fact that this script is called twice,
# once for the whole system and once for the QM subsystem
# Since QM and MM atoms are ordered, this should not be hard

# In principle, here we could also implement link atoms I guess

# prepare_pdb template.pdb > $input.pdb

# Run OPENMM script
python openmm_engrap.py > ../../engrad_mm.dat.$ibead
