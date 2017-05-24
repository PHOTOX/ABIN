#!/usr/bin/env python3

import math
import argparse, sys

# Some usefull string constants
BL = "%BLOCK"
ENDBL = "%ENDBLOCK"
LTCART = "LATTICE_CART"
FRAC = "POSITIONS_FRAC"

def read_cmd():
   """Function for reading command line options."""
   desc = "A simple program for converting a CASTEP cell file to XYZ coordinates."
   parser = argparse.ArgumentParser(description=desc)
   parser.add_argument('-c',dest='cfile',required = True, help='Input cell file.')
   return parser.parse_args()


def vector_norm(vector):
   assert (len(vector)>0),"Empty vector passed!"
   norm = 0.0
   for v in vector:
      norm += v*v
   return math.sqrt(norm)


def dotprod(vec1, vec2):
   assert (len(vec1)==len(vec2)), "Cannot do dotproduct for vectors of different lengths!"
   dot = 0.0
   for i in range(len(vec1)):
      dot += vec1[i]*vec2[i]
   return dot

def read_frac_coords(xfile, atoms, frac):
   vec = [0.0,0.0,0.0]
   read_frac = False
   with open(xfile, "r") as f:
      for l in f:
         ls = l.split()
         if len(ls) < 2:
            continue
         ls0 = ls[0].upper()
         ls1 = ls[1].upper()

         if ls0 == BL and ls1 == FRAC:
            # Start reading fractional coordinates
            read_frac = True
            continue

         if ls0 == ENDBL and ls1 == FRAC:
            # Stop reading fractional coordinates
            read_frac = False
            return

         if read_frac:
            # Sanity check
            if len(ls) < 4:
               print("ERROR during reading of fractional coordinates!")
               print("The offending line is: ", l)
            at, vec[0], vec[1], vec[2] = ls
            atoms.append(at)
            frac.append( (float(vec[0]), float(vec[1]), float(vec[2]) ) )


def read_lattice_cart(cfile, h):
   with open(cfile,"r") as f:
      read_lattice = 4
      for l in f:
         ls = l.split()
         if len(ls) < 2:
            continue
         ls0 = ls[0].upper()
         ls1 = ls[1].upper()

         if ls0 == BL and ls1 == LTCART:
            # Start reading lattice parameters
            read_lattice = 0
            continue
 
         if read_lattice < 3:
            xyz = l.split()
            # A little sanity check
            if len(xyz) < 3:
               print("ERROR when reading lattice parameters from file ", cell_file)
               print("Offending line: ", l)
               sys.exit(1)

            h.append( (float(xyz[0]), float(xyz[1]), float(xyz[2]) ) )
            read_lattice += 1
            continue
 
         # Sanity check
         if read_lattice == 3 and (ls0 != ENDBL or ls1 != LTCART) :
               print("ERROR when reading lattice parameters from file ", cell_file)
               print("Offending line: ", l)
               sys.exit(1)
         else:
            return

def print_xyz_coords(atoms, xyz):
   print( len(xyz) )
   print("")
   i=0
   for h in xyz:
      print('%s   % 18.15f  % 18.15f  % 18.15f' % (atoms[i],h[0],h[1],h[2]))
      i += 1

def convert_frac2xyz(h, frac_coords, xyz):
   for f in frac_coords:
      x = h[0][0]*f[0] + h[1][0]*f[1] + h[2][0]*f[2]
      y = h[0][1]*f[0] + h[1][1]*f[1] + h[2][1]*f[2]
      z = h[0][2]*f[0] + h[1][2]*f[1] + h[2][2]*f[2]
      xyz.append( (x,y,z) )


# MAIN CODE #

opts = read_cmd()
# Input cell file
cell_file = opts.cfile

# Lattice vectors
h = []

#XYZ coordinates
xyz = []

# Atom types
atoms = []

# Fractional coordinates
frac_coords = []

read_lattice_cart(cell_file, h)

read_frac_coords(cell_file, atoms, frac_coords)

convert_frac2xyz(h, frac_coords, xyz)

print_xyz_coords(atoms, xyz)

