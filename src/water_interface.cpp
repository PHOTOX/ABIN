// Copyright 2013 Volodymyr Babin <vb27606@gmail.com>
//
// This is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// The code is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You can find a copy of the GNU General Public License at
// http://www.gnu.org/licenses/.

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../water_potentials/ttm2f.h"
#include "../water_potentials/ttm3f.h"
#include "../water_potentials/ttm4f.h"

#include "../water_potentials/qtip4pf.h"

extern"C" {
#ifdef USE_CP2K
  void force_water(
#else
  void force_water_(
#endif
    const double *x,
    const double *y,
    const double *z,
    double *fx,
    double *fy,
    double *fz,
    double *eclas,
    const int *natom,
    const int* nwalk, 
    const int *watpot)
  {
   const int nwater = *natom / 3;
   double E;
   double grd[9 * nwater];
   double crd[9 * nwater];

   // conversion constants
   const double ANG = 1.889726132873;
   const double AUTOKCAL = 627.50946943;
   const double FAC = 1 / ANG / AUTOKCAL;


   h2o::qtip4pf pot1;
   h2o::ttm2f pot2;
   h2o::ttm3f pot3;
   h2o::ttm4f pot4;

   for (int iw=0;iw < *nwalk;iw++) {

      // Convert to Angstroms
      for (int iat=0; iat < *natom;iat++) {
         crd[3*iat] = x[iw*(*natom)+iat] / ANG;
         crd[3*iat+1] = y[iw*(*natom)+iat] / ANG;
         crd[3*iat+2] = z[iw*(*natom)+iat] / ANG;
      }

      switch ( *watpot) {
      case 1:
         E = pot1(nwater, crd, grd);
         break;
      case 2:
         E = pot2(nwater, crd, grd);
         break;
      case 3:
         E = pot3(nwater, crd, grd);
         break;
      case 4:
         E = pot4(nwater, crd, grd);
         break;
      default:
         std::cerr << "Error: Parameter watpot out of range." << std::endl;
         return;
        break;
      }

      *eclas += E / AUTOKCAL;

      // Convert forces to atomic units (for ABIN)
      for (int iat=0; iat < *natom;iat++) {
         fx[iw*(*natom)+iat]=-grd[3*iat]*FAC;
         fy[iw*(*natom)+iat]=-grd[1+3*iat]*FAC;
         fz[iw*(*natom)+iat]=-grd[2+3*iat]*FAC;
      }

   }

}

// externC
}
