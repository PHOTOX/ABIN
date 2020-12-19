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
#endif /* HAVE_CONFIG_H */

#include <cmath>
#include <cassert>

#include <algorithm>

#include "qtip4pf.h"
#include "constants.h"

//
// this is q-TIP4P/F : http://jcp.aip.org/resource/1/jcpsa6/v131/i2/p024501_s1
//

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

using ttm::CHARGECON;

//----------------------------------------------------------------------------//

const double qM = -1.1128*CHARGECON;
const double qH = -qM/2;

const double rOHeq = 0.9419; // A
const double aHOHeq = 107.4; // degree

const double epsilon = 0.1852; // kcal/mol
const double sigma = 3.1589; // A

const double alpha_r = 2.287; // A-1
const double D_r = 116.09; // kcal/mol
const double k_theta = 87.85; // kcal/mol/rad^2

//----------------------------------------------------------------------------//

const double gammaM = 0.73612;
const double gamma1 = (1.0 - gammaM)/2;

inline void compute_M_site_crd
    (const double* O, const double* H1, const double* H2, double* M)
{
    for (size_t i = 0; i < 3; ++i)
        M[i] = gammaM*O[i] + gamma1*(H1[i] + H2[i]);
}

//----------------------------------------------------------------------------//

double monomer(const double* w, double* grad)
{
    double d1[3], d2[3];
    double r1(0), r2(0), d12(0);

    for (size_t i = 0; i < 3; ++i) {
       d1[i] = w[i + 3] - w[i]; // H1 - O
       r1 += d1[i]*d1[i];

       d2[i] = w[i + 6] - w[i]; // H2 - O
       r2 += d2[i]*d2[i];

       d12 += d1[i]*d2[i];
    }

    r1 = std::sqrt(r1);
    r2 = std::sqrt(r2);

    const double x1 = alpha_r*(r1 - rOHeq);
    const double x2 = alpha_r*(r2 - rOHeq);

    const double V1 = x1*x1*(1.0 - x1*(1.0 - (7.0/12.0)*x1))*D_r;
    const double V2 = x2*x2*(1.0 - x2*(1.0 - (7.0/12.0)*x2))*D_r;

    const double r12 = 1.0/r1/r2;
    const double xi = d12*r12;
    const double theta = std::acos(xi);
    const double d_theta = theta - aHOHeq*M_PI/180.0;

    if (grad) {
        const double dV1 = alpha_r*D_r*x1*(2.0 - x1*(3.0 - (7.0/3.0)*x1))/r1;
        const double dV2 = alpha_r*D_r*x2*(2.0 - x2*(3.0 - (7.0/3.0)*x2))/r2;

        const double factor = - k_theta*d_theta/std::sqrt(1.0 - xi*xi);

        const double t0 = factor*r12;
        const double t1 = dV1 - factor*xi/r1/r1;
        const double t2 = dV2 - factor*xi/r2/r2;

        for (size_t i = 0; i < 3; ++i) {
            const double dA1 = t0*d2[i] + t1*d1[i];
            const double dA2 = t0*d1[i] + t2*d2[i];

            grad[i] = -(dA1 + dA2);
            grad[i + 3] = dA1;
            grad[i + 6] = dA2;
        }
    }

    return V1 + V2 + 0.5*k_theta*d_theta*d_theta;
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o {

//----------------------------------------------------------------------------//

qtip4pf::qtip4pf()
{
    m_nw = 0;
    m_mem = 0;
}

//----------------------------------------------------------------------------//

qtip4pf::~qtip4pf()
{
    delete[] m_mem;
}

//----------------------------------------------------------------------------//

void qtip4pf::allocate(size_t nw)
{
    if (m_nw > nw || nw == 0)
        return;

    delete[] m_mem;

    // see operator() below
    m_mem = new double[3*nw + 3*nw];

    m_nw = nw;
}

//----------------------------------------------------------------------------//

double qtip4pf::operator()(size_t nw, const double* crd)
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers
    double* msite = m_mem;

    // compute M-sites and 1-body terms
    double Eint(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        compute_M_site_crd(crd + i9, crd + i9 + 3, crd + i9 + 6, msite + i3);
        Eint += monomer(crd + i9, 0);
    }

    // vdW and electrostatics
    double Evdw(0), Eelec(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        // Hydrogens - M-sites
        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            const size_t j3 = 3*j;

            for (size_t l = 1; l < 3; ++l) {
                const size_t j3l = j3 + l;
                const size_t jh = 3*j3l;

                double Rsq(0);
                for (size_t k = 0; k < 3; ++k) {
                    const double dx = msite[i3 + k] - crd[jh + k];
                    Rsq += dx*dx;
                }

                Eelec += qH*qM/std::sqrt(Rsq);
            }
        }

        // M-sites, Oxygens and Hydrogen - Hydrogen

        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j3 = 3*j;
            const size_t j9 = 3*j3;

            // H-H
            for (size_t a = 1; a < 3; ++a) {
                const size_t i3a = i3 + a;

                for (size_t b = 1; b < 3; ++b) {
                    const size_t j3b = j3 + b;

                    double Rsq(0);
                    for (size_t k = 0; k < 3; ++k) {
                        const double dx = crd[3*i3a + k] - crd[3*j3b + k];
                        Rsq += dx*dx;
                    }

                    Eelec += qH*qH/std::sqrt(Rsq);
                }
            }

            // O-O distance
            double Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                const double dx = crd[i9 + k] - crd[j9 + k];
                Rsq += dx*dx;
            }

            Rsq = sigma*sigma/Rsq;

            const double dR6 = Rsq*Rsq*Rsq;

            Evdw += 4*epsilon*dR6*(dR6 - 1.0);

            // M-M distance
            Rsq = 0.0;
            for (size_t k = 0; k < 3; ++k) {
                const double dx = msite[i3 + k] - msite[j3 + k];
                Rsq += dx*dx;
            }

            Eelec += qM*qM/std::sqrt(Rsq);
        }
    }

    return Eint + Evdw + Eelec;
}

//----------------------------------------------------------------------------//

double qtip4pf::operator()(size_t nw, const double* crd, double* grad)
{
    if (grad == 0)
        return operator()(nw, crd);

    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers
    double* msite = m_mem;
    double* dM = msite + 3*nw;

    // compute M-sites and 1-body terms
    double Eint(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        compute_M_site_crd(crd + i9, crd + i9 + 3, crd + i9 + 6, msite + i3);
        Eint += monomer(crd + i9, grad + i9);
    }

    std::fill(dM, dM + 3*nw, 0.0);

    // vdW and electrostatics
    double Evdw(0), Eelec(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        // Hydrogens - M-sites
        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            const size_t j3 = 3*j;

            for (size_t l = 1; l < 3; ++l) {
                const size_t j3l = j3 + l;
                const size_t jh = 3*j3l;

                double dR[3], Rsq(0);
                for (size_t k = 0; k < 3; ++k) {
                    dR[k] = msite[i3 + k] - crd[jh + k];
                    Rsq += dR[k]*dR[k];
                }

                const double r1 = 1.0/std::sqrt(Rsq);
                const double E1 = qH*qM*r1;
                const double g1 = E1*r1*r1;

                Eelec += E1;

                for (size_t k = 0; k < 3; ++k) {
                    dM[i3 + k] -= g1*dR[k];
                    grad[jh + k] += g1*dR[k];
                }
            }
        }

        // M-sites, Oxygens and Hydrogen - Hydrogen

        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j3 = 3*j;
            const size_t j9 = 3*j3;

            // H-H
            for (size_t a = 1; a < 3; ++a) {
                const size_t i3a = i3 + a;

                for (size_t b = 1; b < 3; ++b) {
                    const size_t j3b = j3 + b;

                    double dR[3], Rsq(0);
                    for (size_t k = 0; k < 3; ++k) {
                        dR[k] = crd[3*i3a + k] - crd[3*j3b + k];
                        Rsq += dR[k]*dR[k];
                    }

                    const double r1 = 1.0/std::sqrt(Rsq);
                    const double E1 = qH*qH*r1;
                    const double g1 = E1*r1*r1;

                    Eelec += E1;

                    for (size_t k = 0; k < 3; ++k) {
                        grad[3*i3a + k] -= g1*dR[k];
                        grad[3*j3b + k] += g1*dR[k];
                    }
                }
            }

            // O-O distance
            double dR[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                dR[k] = crd[i9 + k] - crd[j9 + k];
                Rsq += dR[k]*dR[k];
            }

            const double sRsq = sigma*sigma/Rsq;
            const double dR6 = sRsq*sRsq*sRsq;

            Evdw += 4*epsilon*dR6*(dR6 - 1.0);
            const double g_vdw = 24*epsilon*dR6*(2*dR6 - 1.0)/Rsq;

            for (size_t k = 0; k < 3; ++k) {
                grad[i9 + k] -= g_vdw*dR[k];
                grad[j9 + k] += g_vdw*dR[k];
            }

            // M-M distance
            Rsq = 0.0;
            for (size_t k = 0; k < 3; ++k) {
                dR[k] = msite[i3 + k] - msite[j3 + k];
                Rsq += dR[k]*dR[k];
            }

            const double r1 = 1.0/std::sqrt(Rsq);
            const double E1 = qM*qM*r1;
            const double g1 = E1*r1*r1;

            Eelec += E1;

            for (size_t k = 0; k < 3; ++k) {
                dM[i3 + k] -= g1*dR[k];
                dM[j3 + k] += g1*dR[k];
            }
        }
    }

    // distribute M-site gradients

    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        for (size_t k = 0; k < 3; ++k) {
            grad[i9 + 0 + k] += gammaM*dM[i3 + k]; // O
            grad[i9 + 3 + k] += gamma1*dM[i3 + k]; // H
            grad[i9 + 6 + k] += gamma1*dM[i3 + k]; // H
        }
    }

    return Eint + Evdw + Eelec;
}

//----------------------------------------------------------------------------//

} // namespace h2o

////////////////////////////////////////////////////////////////////////////////
