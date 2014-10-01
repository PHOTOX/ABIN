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

#include "ttm2f.h"
#include "ttm3f.h"
#include "ttm4f.h"

#include "qtip4pf.h"

////////////////////////////////////////////////////////////////////////////////
int force_water(const double x[],const double  y[],const double z[], double fx[],double fy[],double fz[], \
      double E,const int nw, const int watpot)
{
   double grd[9*nw];

   double crd[9*nw];

   h2o::qtip4pf pot1;
   h2o::ttm2f pot2;
   h2o::ttm3f pot3;
   h2o::ttm4f pot4;

   switch (watpot) {
   case 1:
      E = pot1(nw, crd, grd);
      break;
   case 2:
      E = pot2(nw, crd, grd);
      break;
   case 3:
      E = pot3(nw, crd, grd);
      break;
   case 4:
      E = pot4(nw, crd, grd);
      break;
   default:
      std::cerr << "Error: Parameter watpot out of range." << std::endl;
      return 1;
      break;
   }


   return 0;

}

////////////////////////////////////////////////////////////////////////////////
