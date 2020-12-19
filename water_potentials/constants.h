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

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace ttm {

// 1D \approx 0.20822678 e*A
const double DEBYE = 1.0/2.081943416923e-01;

// my $Fr = 4.8032042712e-10; # 1e in ESU
// my $cal_to_J = 4.184; # 1 cal in J
// my $Na = 6.0221417930e+23;
//
// my $R = 1.0e-8; # 1A in cm
//
// my $E_cc = $Fr*$Fr/$R; # in ergs; 1 erg = 10e-7 J
//
// $E_cc *= 1.0e-7; # now in J
// $E_cc /= $cal_to_J; # now in cal
// $E_cc *= $Na; # now in cal/mol
// $E_cc /= 1000; # now in kcal/mol
//
// print "CHARGECON = ", sqrt($E_cc), "\n";

const double CHARGECON = 1.822261544733e+01;

} // namespace ttm

#endif // CONSTANTS_H
