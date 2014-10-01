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

#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <algorithm>

#include "ps.h"
#include "smear.h"
#include "constants.h"

#include "ttm4f.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

// vdw
const double A6  = -0.126503e+5;
const double A8  =  0.526347e+6;
const double A10 = -0.964270e+7;
const double A12 =  0.877792e+8;
const double A14 = -0.367476e+9;
const double A16 =  0.572395e+9;

// M-site positioning
const double gammaM = 0.426706882;
const double gamma1 = 1.0 - gammaM;
const double gamma2 = gammaM/2;

// polarizability
const double polarO = 1.310;
const double polarH = 0.294;
const double polarM = 0.0;

// Thole damping factors
const double polfacO = polarO;
const double polfacH = polarH;
const double polfacM = polarO;

using ttm::CHARGECON;

//----------------------------------------------------------------------------//

inline void compute_M_site_crd
    (const double O[3], const double H1[3], const double H2[3], double M[3])
{
    for (size_t i = 0; i < 3; ++i)
        M[i] = gamma1*O[i] + gamma2*(H1[i] + H2[i]);
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o {

//----------------------------------------------------------------------------//

ttm4f::ttm4f()
{
    m_nw = 0;
    m_memory = 0;
}

//----------------------------------------------------------------------------//

ttm4f::~ttm4f()
{
    delete[] m_memory;
}

//----------------------------------------------------------------------------//

void ttm4f::allocate(size_t nw)
{
    if (m_nw > nw || nw == 0)
        return;

    delete[] m_memory;

    m_nw = nw;
    m_memory = new double[4*m_nw // charge
                        + 4*m_nw  // polfac
                        + 4*m_nw   // polar
                        + 3*4*m_nw  // xyz
                        + 3*3*3*m_nw // grdq
                        + 3*4*m_nw];  // dE

    // setup exclusions
    m_excluded.clear();
    for (size_t i = 0; i < m_nw; ++i) {
        const size_t i4 = 4*i;
        m_excluded.insert(std::make_pair(i4 + 0, i4 + 1)); // O - H1
        m_excluded.insert(std::make_pair(i4 + 0, i4 + 2)); // O - H2
        m_excluded.insert(std::make_pair(i4 + 0, i4 + 3)); // O - M
        m_excluded.insert(std::make_pair(i4 + 1, i4 + 2)); // H1 - H2
        m_excluded.insert(std::make_pair(i4 + 1, i4 + 3)); // H1 - M
        m_excluded.insert(std::make_pair(i4 + 2, i4 + 3)); // H2 - M
    }
}

//----------------------------------------------------------------------------//

double ttm4f::operator()(size_t nw, const double* crd, double* grad)
{
    assert(crd != 0 || nw == 0);

    allocate(nw);

    // setup pointers

    double* charge = m_memory;
    double* polfac = charge + 4*nw;
    double* polar = polfac + 4*nw;
    double* xyz = polar + 4*nw;
    double* grdq = xyz + 3*4*nw;
    double* dE = grdq + 3*3*3*nw;

    // populate polfac and polar
    for (size_t n = 0; n < nw; ++n) {
        const size_t n4 = 4*n;

        polfac[n4 + 0] = polfacO;
        polfac[n4 + 1] = polfacH;
        polfac[n4 + 2] = polfacH;
        polfac[n4 + 3] = polfacM;

        polar[n4 + 0] = polarO;
        polar[n4 + 1] = polarH;
        polar[n4 + 2] = polarH;
        polar[n4 + 3] = polarM;
    }


    // calculate coordinates of the M-sites and
    // pack the coordinates in xyz as O H H M O H H M ...

    for (size_t i = 0; i < nw; ++i) {
        // O H H
        std::copy(crd + 3*3*i, crd + 3*3*i + 9, xyz + 4*3*i);

        // and M
        compute_M_site_crd(crd + 3*3*i, crd + 3*(3*i + 1), crd + 3*(3*i + 2),
                           xyz + 4*3*i + 9);
    }

    // calculate the INTRA-molecular energy and Dipole Moment
    // Surface (Partridge-Schwenke)

    double Eint = 0.0;
    for (size_t n = 0; n < nw; ++n) {
        double q3[3], dq3[27];
        const size_t n9 = 9*n;

        if (grad) {
            Eint += ps::pot_nasa(crd + n9, grad + n9);
            ps::dms_nasa(0.0, 0.0, 0.0,
                         crd + n9, q3, dq3, false);
        } else {
            Eint += ps::pot_nasa(crd + n9, 0);
            ps::dms_nasa(0.0, 0.0, 0.0,
                         crd + n9, q3, 0, false);
        }

        // TTM2.1-F assignment
        const double tmp = 0.5*gammaM/(1.0 - gammaM);

        charge[4*n + 0] = 0.0;                        // O
        charge[4*n + 1] = q3[1] + tmp*(q3[1] + q3[2]); // H1
        charge[4*n + 2] = q3[2] + tmp*(q3[1] + q3[2]); // H2
        charge[4*n + 3] = q3[0]/(1.0 - gammaM);       // M

        if (grad == 0)
            continue;

#define DQ3(i,j,k) dq3[k + 3*(j + 3*i)]
#define GRDQ(i,j,k) grdq[k + 3*(j + 3*(i + 3*n))]

        for (size_t k = 0; k < 3; ++k) {
            GRDQ(0, 0, k) = DQ3(0, 0, k) + tmp*(DQ3(0, 0, k) + DQ3(0, 1, k));
            GRDQ(1, 0, k) = DQ3(1, 0, k) + tmp*(DQ3(1, 0, k) + DQ3(1, 1, k));
            GRDQ(2, 0, k) = DQ3(2, 0, k) + tmp*(DQ3(2, 0, k) + DQ3(2, 1, k));

            GRDQ(0, 1, k) = DQ3(0, 1, k) + tmp*(DQ3(0, 1, k) + DQ3(0, 0, k));
            GRDQ(1, 1, k) = DQ3(1, 1, k) + tmp*(DQ3(1, 1, k) + DQ3(1, 0, k));
            GRDQ(2, 1, k) = DQ3(2, 1, k) + tmp*(DQ3(2, 1, k) + DQ3(2, 0, k));

            GRDQ(0, 2, k) = DQ3(0, 2, k) - 2*tmp*(DQ3(0, 0, k) + DQ3(0, 1, k));
            GRDQ(1, 2, k) = DQ3(1, 2, k) - 2*tmp*(DQ3(1, 0, k) + DQ3(1, 1, k));
            GRDQ(2, 2, k) = DQ3(2, 2, k) - 2*tmp*(DQ3(2, 0, k) + DQ3(2, 1, k));
        }
    }

    for (size_t n = 0; n < 4*nw; ++n)
        charge[n] *= CHARGECON;

    if (grad) {
        for (size_t n = 0; n < 27*nw; ++n)
            grdq[n] *= CHARGECON;
    }

    //---------------------------------------------------------------!
    // Calculate the vdW interactions for all atoms                  !
    //---------------------------------------------------------------!

    double Evdw(0);

    for (size_t i = 0; i < nw - 1; ++i) {
        for (size_t j = i + 1; j < nw; ++j) {

            // Calculate Oxygen-Oxygen interactions
            double dRijsq = 0.0, Rij[3];
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = crd[3*3*i + k] - crd[3*3*j + k];
                dRijsq += Rij[k]*Rij[k];
            }

            // Add Evdw van der Waal energy for O-O interaction
            const double x = 1/dRijsq;
            Evdw += x*x*x*(A6 + x*(A8 
                              + x*(A10 + x*(A12 + x*(A14 + x*A16)))));

            if (grad) {
                const double tmp =
                    x*x*x*x*(6*A6 + x*(8*A8 + x*(10*A10
                             + x*(12*A12 + x*(14*A14 + x*16*A16)))));

                for (size_t k = 0; k < 3; ++k) {
                    grad[3*3*i + k] -= tmp*Rij[k];
                    grad[3*3*j + k] += tmp*Rij[k];
                }
            } // grad
        }
    }

    ttm::smear_ttm4 smr;

    if (grad == 0)
        return Eint + Evdw
          + m_electrostatics(4*nw, charge, polfac,
                             polar, xyz, m_excluded, smr, 0);

    //-------------------------------------------------------------------------!
    // Calculate the remaining part of the derivatives                         !
    //-------------------------------------------------------------------------!

    const double Eelec =
        m_electrostatics(4*nw, charge, polfac,
                         polar, xyz, m_excluded, smr, dE);

    assert(m_electrostatics.dipoles_converged());

    // add electrostatic derivatives
    for (size_t n = 0; n < nw; ++n) {
        // O H H
        for (size_t k = 0; k < 9; ++k)
            grad[9*n + k] += dE[12*n + k];

        // redistribute the M-site derivatives
        for (size_t k = 0; k < 3; ++k) {
            grad[9*n + 0 + k] += (1.0 - gammaM)*dE[12*n + 9 + k]; // O
            grad[9*n + 3 + k] +=     0.5*gammaM*dE[12*n + 9 + k]; // H
            grad[9*n + 6 + k] +=     0.5*gammaM*dE[12*n + 9 + k]; // H
        }
    }

    const double* phi = m_electrostatics.phi();

    // derivatives from the adjustable charges of the NASA's PES
    for (size_t n = 0; n < nw; ++n) {
        const size_t io  = 9*n + 0;
        const size_t ih1 = 9*n + 3;
        const size_t ih2 = 9*n + 6;

        for (size_t k = 0; k < 3; ++k) {
            grad[ih1 + k] += GRDQ(0, 0, k)*phi[4*n + 1]  // phi(h1)
                           + GRDQ(0, 1, k)*phi[4*n + 2]  // phi(h2)
                           + GRDQ(0, 2, k)*phi[4*n + 3]; // phi(M)

            grad[ih2 + k] += GRDQ(1, 0, k)*phi[4*n + 1]  // phi(h1)
                           + GRDQ(1, 1, k)*phi[4*n + 2]  // phi(h2)
                           + GRDQ(1, 2, k)*phi[4*n + 3]; // phi(M)

            grad[io + k] += GRDQ(2, 0, k)*phi[4*n + 1]  // phi(h1)
                          + GRDQ(2, 1, k)*phi[4*n + 2]  // phi(h2)
                          + GRDQ(2, 2, k)*phi[4*n + 3]; // phi(M)
        }
    }

    return Eint + Evdw + Eelec;
}

//----------------------------------------------------------------------------//

} // namespace h2o

////////////////////////////////////////////////////////////////////////////////
