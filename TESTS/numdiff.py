#!/usr/bin/env python
from optparse import OptionParser
import decimal
from sys import argv, exit 
# Workaround for STDIN problem for CP2K tests
import readline

# This scripts is compares differences between files
# while disregarding insignificant numerical differences.

# This script parses output from command:
# diff -y file1 file2

# It returns 1 when difference is found and 0 otherwise

decimal.getcontext().prec = 17

# TODO: Relative delta
DELTA=1e-16
#decimal.getcontext().rounding = "ROUND_DOWN"

usage = "USAGE: %prog input_diff_file"
parser = OptionParser(usage)
options, args = parser.parse_args(argv[1:])
try:
   inpfile = args[0]
except:
   print("Error: You need to provide input file as an argument.")
   print("Type -h for help")
   exit(1)

print("Comparing numerical differences in file "+inpfile)

def failed(reference, test):
   delta = float(reference) - float(test)
   # Let's just not report small differences at all,
   # It just creates a noise
   # TODO: Testing only 2 steps should allow to tighten this up significantly
   if delta < DELTA:
      return
   print('>>>>')
   print('Significant numerical difference in file %s' % inpfile)
   print('Reference: ' + reference)
   print('Test: ' + test)
   print('Delta = ', delta)
   exit(1)

def compare(numbers1, numbers2):
    """Compare each element of numbers1 and numbers2"""
    for i in range(len(numbers1)):
       if numbers1[i] == numbers2[i]:
          continue

       try:
          # Rounding to significant precision
          dec1 = +decimal.Decimal(numbers1[i])
          dec2 = +decimal.Decimal(numbers2[i]) #*decimal.Decimal(1)
       except decimal.InvalidOperation:
          print("Non-numerical difference found!")
          print("line1 = ", numbers1[i])
          print("line2 = ", numbers2[i])
          exit(1)

       if dec1 != dec2:
          # Compare sign, exponents, and digits separately
          sign1, digits1, exp1 = dec1.as_tuple() 
          sign2, digits2, exp2 = dec2.as_tuple() 
          if sign1 != sign2 or exp1 != exp2:
             failed(numbers1[i], numbers2[i])
             continue

          ints1=""
          ints2=""
          # Compare last 10 significant digits
          # TODO: Why aren't we comparing all digits?
          for j in digits1[-10:]:
             ints1 += str(j)
          for j in digits2[-10:]:
             ints2 += str(j)
          ints1 = int(ints1)
          ints2 = int(ints2)
          # This allows to get rid off insignificant rounding errors
          if abs(ints1 - ints2) != 1:
             failed(numbers1[i], numbers2[i])


with open(inpfile, 'r') as f:
   # If the file is empty something is wrong 
   # (e.g. file for comparison was not even generated)
   if len(f.read()) == 0:
      print("Blank file encountered.")
      exit(1)
   f.seek(0)

   for line in f:
      split = line.split()
      # Exit if there are non-numerical differences
      if "<" in split or ">" in split:
           print("There's an extra line!")
           print(line)
           exit(1)
      # Skip line if there is no difference
      if '|' not in split:
          continue

      split = line.split('|')
      diff1 = split[0].split()
      diff2 = split[1].split()
      if len(diff1) != len(diff2):
          print("Number of columns differ!")
          exit(1)
      compare(diff1, diff2)

exit(0)

