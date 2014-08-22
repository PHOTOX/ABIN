#!/usr/bin/awk -f
# Tiny awk script for computing average and std. deviation from one column of data,
# It needs variable "column" defined by the caller
# Should be called as:
# prum.awk -v column=$column inputfile
BEGIN{x=0.0;xx=0.0}
{
 x=x+$column
 xx=xx+$column*$column
}
END{
# print (x/NR,xx/NR,xx/NR-(x/NR)^2,sqrt(xx/NR-(x/NR)^2))
 print "# Column   Average    Std. deviation"
 print "#",column, x/NR,sqrt(xx/NR-(x/NR)^2)
}

