#!/bin/bash
#
if [ -z "$1" ];then
	echo "USAGE: $0 [molpro.output]"
	echo 
	echo "Script for extracting geometry from molpro output files for visualization in molden."
	echo "Molden data are stored in temporary file ""temp"". Quit molden with ""Ctrl+C"" in order to access this file."
	exit 1
fi

if [ ! -e $1 ];then
	echo "File $1 does not exists."
	echo "USAGE: $0 [molpro.output]"
	exit 1
fi
#reading XYZ geometry in molpro output,conversion from bohr to angstroms
awk '
{
	if ( ($1 == "NR" && $2 == "ATOM")||($1 == "Nr" && $2 == "Atom") ) {
		n = 0
		getline
		getline
		while ( $1 >= 0 ) {
			n++
			a[n] = $2 
			x[n] = $4/1.8897 
			y[n] = $5/1.8897 
			z[n] = $6/1.8897
			getline
		}
		print "  "n
		print " "
		for (k=1;k<=n;k++)
		       printf "%2s   % 4.8f  % 4.8f % 4.8f\n",a[k],x[k],y[k],z[k]	
	}

}' $1 
#################### NOW MOLDEN DOES ITS JOB
if [ $? != "0" ];then
        echo "An error encountered.Exiting..."
	exit $? 
fi

