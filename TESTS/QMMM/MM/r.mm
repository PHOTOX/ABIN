#!/bin/bash
cd MM
ORCAEXE=orca


timestep=$1
ibead=$2
input=input$ibead
geom=../geom_mm.dat.$ibead

natom=`cat $geom | wc -l`

rm -f $input.engrad $input.gbw

cat > $input << EOF
# in the following line, specify basis sets and method,
# basis sets for RI approximations are basis/J
# PAL2
# ! RI BP86 SVP SVP/J
! PM3
! ENGRAD  #TightSCF
* xyz 0 1
EOF
########END OF USE MODIFICATIONS###################

cat $geom >> $input
echo '*' >>$input

$ORCAEXE $input &> $input.out
################################
cp $input.out $input.out.old

### EXTRACTING ENERGY AND FORCES
awk -v natom="$natom" '{if ($5=="energy") {getline;getline;print $1}
if ($4=="gradient") {getline;
	for(i=1;i<=natom;i++) {
		getline; x=$1;getline;y=$1;getline;print x,y,$1
 }}}' $input.engrad > ../engrad_mm.dat.$ibead



