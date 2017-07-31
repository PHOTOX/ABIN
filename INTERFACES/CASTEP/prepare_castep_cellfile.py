#!/usr/bin/env python
from __future__ import print_function
import math
import numpy as np
import argparse, sys

# Some usefull string constants
BL = "%BLOCK"
ENDBL = "%ENDBLOCK"
LTCART = "LATTICE_CART"
FRAC = "POSITIONS_FRAC"


def read_cmd():
   """Function for reading command line options."""
   desc = "A simple program for converting XYZ coordinates to a CASTEP cell file."
   parser = argparse.ArgumentParser(description=desc)
   parser.add_argument('-c',dest='cfile',required = True, help='Input template cell file.')
   parser.add_argument('-x',dest='xfile',required = True, help='Input XYZ file.')
   return parser.parse_args()


def dotprod(vec1, vec2):
   assert (len(vec1)==len(vec2)), "Cannot do dotproduct for vectors of different lengths!"
   dot = 0.0
   for i in range(len(vec1)):
      dot += vec1[i]*vec2[i]
   return dot


def print_lattice_cart(h_matrix):
   print(BL, LTCART)
   for h in h_matrix:
      print('% 18.15f  % 18.15f  % 18.15f' % (h[0],h[1],h[2]))
   print(ENDBL, LTCART)
   print("")


def print_frac_coords(atoms, fcoords):
   print(BL, FRAC)
   i=0
   for h in fcoords:
      print('%s   % 18.15f  % 18.15f  % 18.15f' % (atoms[i],h[0],h[1],h[2]))
      i += 1

   print(ENDBL, FRAC)
   print("")


def print_cellfile(h, frac, atoms,  cfile):
   print_lattice_cart(h)
   print_frac_coords(atoms, frac)

   # Now print the rest of the template cell file
   with open(cfile,"r") as f:
      print_line = True

      for l in f:
         ls = l.split()
         if len(ls) == 0:
            continue
         else:
            ls0 = ls[0].upper()

         if len(ls) > 1:
            ls1 = ls[1].upper()
         else:
            ls1 = ""

         # Print all blocks except for lattice cart and coordinates
         if ls0==BL and (ls1 == LTCART or ls1 == FRAC):
            print_line = False
            continue

         if not print_line and ls0==ENDBL:
            print_line = True
            continue

         if print_line:
            print(l.split("\n")[0])


def read_xyz_coords(xfile, atoms, xyz):
   vec = [0.0,0.0,0.0]
   with open(xfile, "r") as f:
      for l in f:
         at, vec[0], vec[1], vec[2] = l.split()
         atoms.append(at)
         xyz.append( (float(vec[0]), float(vec[1]), float(vec[2]) ) )


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


def convert_xyz2frac(h, xyz, frac_coords):
   sigma_a = np.cross(h[1],h[2])
   sigma_b = np.cross(h[2],h[0])
   sigma_c = np.cross(h[0],h[1])
   Omega = dotprod(h[0],sigma_a)
   
   cart2frac = [sigma_a/Omega, sigma_b/Omega, sigma_c/Omega]
   
   cart = [0.0,0.0,0.0]
   for i in range(len(xyz)):
            cart = [float(xyz[i][0]), float(xyz[i][1]), float(xyz[i][2])]
            frac = [0.0,0.0,0.0]
            for i in range(0,3):
                for j in range(0,3):
                    frac[i] = cart2frac[i][j]*cart[j] + frac [i]
            frac_coords.append((frac[0], frac[1], frac[2]))

# MAIN CODE #

opts = read_cmd()
# Input cell file
cell_file = opts.cfile
# Input XYZ file
xyz_file = opts.xfile

# Lattice vectors
h = []

#XYZ coordinates
xyz = []

# Atom types
atoms = []

# Fractional coordinates
frac_coords = []

read_lattice_cart(cell_file, h)

read_xyz_coords(xyz_file, atoms, xyz)

convert_xyz2frac(h, xyz, frac_coords)

print_cellfile(h, frac_coords, atoms, cell_file)

