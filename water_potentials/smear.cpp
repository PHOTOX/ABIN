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

#include <limits>

#include "gammq.h"
#include "smear.h"

namespace {

const double EPS = std::numeric_limits<double>::epsilon();

} // namespace

namespace ttm {

////////////////////////////////////////////////////////////////////////////////

void smear_ttm3::smear01
    (const double& drsq, const double& pol12, const double& a,
     double &ts0, double &ts1) const
{
    const double dd = std::sqrt(drsq);
    const double dri = 1.0/dd;
    const double drsqi = dri*dri;

    if (pol12 > EPS) {
        const double g23 = std::exp(gammln(2.0/3.0));

        const double AA = std::pow(pol12, 1.0/6.0);
        const double rA = dd/AA;
        const double rA3 = std::pow(rA, 3);
        const double exp1 = std::exp(-a*rA3);
        const double a_cubert = std::pow(a, 1.0/3.0);

        ts0 = (1.0 - exp1 + a_cubert*rA*g23*ttm::gammq(2.0/3.0, a*rA3))*dri;
        ts1 = (1.0 - exp1)*dri*drsqi;
    } else {
        ts0 = dri;
        ts1 = dri*drsqi;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smear_ttm3::smear2
    (const double& drsq, const double& pol12, const double& a,
     double &ts1, double &ts2) const
{
    const double dd = std::sqrt(drsq);
    const double dri = 1.0/dd;
    const double drsqi = dri*dri;

    if (pol12 > EPS) {
        const double AA = std::pow(pol12, 1.0/6.0);
        const double rA = dd/AA;
        const double rA3 = std::pow(rA, 3);
        const double exp1 = std::exp(-a*rA3);

        ts1 = (1.0 - exp1)*dri*drsqi;
        ts2 = (ts1 - exp1*a/std::pow(AA, 3))*drsqi;
    } else {
        ts1 = dri*drsqi;
        ts2 = ts1*drsqi;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smear_ttm3::smear3
    (const double& drsq, const double& pol12, const double& a,
     double &ts1, double &ts2, double &ts3) const
{
    const double dd = std::sqrt(drsq);
    const double dri = 1.0/dd;
    const double drsqi = dri*dri;

    if (pol12 > EPS) {
        const double AA = std::pow(pol12, 1.0/6.0);
        const double rA = dd/AA;
        const double rA3 = std::pow(rA, 3);
        const double exp1 = std::exp(-a*rA3);

        ts1 = (1.0 - exp1)*dri*drsqi;
        ts2 = (ts1 - exp1*a/std::pow(AA, 3))*drsqi;
        ts3 = (ts2 - 0.6*exp1*dd*a*a/pol12)*drsqi;
    } else {
        ts1 = dri*drsqi;
        ts2 = ts1*drsqi;
        ts3 = ts2*drsqi;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smear_ttm4::smear01
    (const double& drsq, const double& pol12, const double& a,
     double &ts0, double &ts1) const
{
    const double dd = std::sqrt(drsq);
    const double dri = 1.0/dd;
    const double drsqi = dri*dri;

    if (pol12 > EPS) {
        const double g34 = std::exp(gammln(3.0/4.0));

        const double AA = std::pow(pol12, 1.0/6.0);
        const double rA = dd/AA;
        const double rA4 = std::pow(rA, 4);
        const double exp1 = std::exp(-a*rA4);
        const double a_mrt = std::pow(a, 1.0/4.0);

        ts0 = (1.0 - exp1 + a_mrt*rA*g34*ttm::gammq(3.0/4.0, a*rA4))*dri;
        ts1 = (1.0 - exp1)*dri*drsqi;
    } else {
        ts0 = dri;
        ts1 = dri*drsqi;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smear_ttm4::smear2
    (const double& drsq, const double& pol12, const double& a,
     double &ts1, double &ts2) const
{
    const double dd = std::sqrt(drsq);
    const double dri = 1.0/dd;
    const double drsqi = dri*dri;

    if (pol12 > EPS) {
        const double AA = std::pow(pol12, 1.0/6.0);
        const double rA = dd/AA;
        const double rA4 = std::pow(rA, 4);
        const double exp1 = std::exp(-a*rA4);

        ts1 = (1.0 - exp1)*dri*drsqi;
        ts2 = (ts1 - (4.0/3.0)*a*exp1*rA4*dri*drsqi)*drsqi;
    } else {
        ts1 = dri*drsqi;
        ts2 = ts1*drsqi;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smear_ttm4::smear3
    (const double& drsq, const double& pol12, const double& a,
     double &ts1, double &ts2, double &ts3) const
{
    const double dd = std::sqrt(drsq);
    const double dri = 1.0/dd;
    const double drsqi = dri*dri;

    if (pol12 > EPS) {
        const double AA = std::pow(pol12, 1.0/6.0);
        const double AA4 = std::pow(AA, 4);
        const double rA = dd/AA;
        const double rA4 = std::pow(rA, 4);
        const double exp1 = std::exp(-a*rA4);

        ts1 = (1.0 - exp1)*dri*drsqi;
        ts2 = (ts1 - (4.0/3.0)*a*exp1*rA4*dri*drsqi)*drsqi;
        ts3 = (ts2 - (4.0/15.0)*dri*a*(4*a*rA4 - 1.0)*exp1/AA4)*drsqi;
    } else {
        ts1 = dri*drsqi;
        ts2 = ts1*drsqi;
        ts3 = ts2*drsqi;
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace ttm
