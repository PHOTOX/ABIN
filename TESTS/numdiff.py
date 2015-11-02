#!/usr/bin/env python
import argparse
from sys import exit
from decimal import *

getcontext().prec = 14

parser = argparse.ArgumentParser()
parser.add_argument('inpfile', metavar='input_file.diff', help='Text file generetad by diff -y command' )
opts = parser.parse_args()
inpfile = opts.inpfile

def compare(numbers1, numbers2):
    """Compare each element of numbers1 and numbers2"""
    for i in range(len(numbers1)):
        dec1 = Decimal(numbers1[i])*Decimal(1)
        dec2 = Decimal(numbers2[i])*Decimal(1)
        if dec1 != dec2:
            print 'Comparison failed'
            print numbers1[i], numbers2[i]
            exit(1)


with open(inpfile, 'r') as f:
   # If the file is empty something is wrong 
   # (e.g. file for comparison was not even generated)
   if len(f.read())==0:
      print("Blank file encoutered.")
      exit(1)
   f.seek(0)
   for line in f:
      # Exit if there are non-numerical differences
      if "<" in line.split() or ">" in line.split():
           print("There's an extra line!")
           print(line)
           exit(1)
      # Skip line if there is no difference
      if '|' not in line.split():
          continue 

      tmp1 = line.split('|')[0]
      tmp2 = line.split('|')[1]
      diff1 = tmp1.split()
      diff2 = tmp2.split()
      if len(diff1) != len(diff2):
          print("Number of columns differ!")
          exit(1)
      compare(diff1, diff2)

exit(0)


