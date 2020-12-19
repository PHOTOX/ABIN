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

#include "io-xyz.h"
#include "xyz-water-utils.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

template <typename POT_TYPE>
double check_gradient(POT_TYPE& pot, size_t nw, double* crd)
{
    const double EPS = 1.0e-5;

    double grd[9*nw], diff(0);

    (void) pot(nw, crd, grd);

    for (size_t i = 0; i < 9*nw; ++i) {
        const double orig = crd[i];

        crd[i] = orig + EPS;
        const double fp = pot(nw, crd);

        crd[i] = orig - EPS;
        const double fm = pot(nw, crd);

        const double grd_fd = 0.5*(fp - fm)/EPS;
        crd[i] = orig;

        diff += std::pow(grd[i] - grd_fd, 2);
    }

    return std::sqrt(diff);
}

//----------------------------------------------------------------------------//

void hline(std::ostream& os)
{
    os << '\n';
    for (size_t n = 0; n < 80; ++n)
        os << '=';
    os << '\n';
}

//----------------------------------------------------------------------------//

template <typename POT_TYPE>
void invoke_potential(POT_TYPE& pot, size_t nw, double* crd)
{
    const double E = pot(nw, crd);
    const double grad_err = check_gradient(pot, nw, crd);

    hline(std::cout);
    std::cout << std::setw(9) << pot.name() << " : E = " << E << " kcal/mol\n"
                 "|analytical_gradient - finite_differences| = "
              << grad_err << " kcal/mol/A\n";
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: test-ttm water.xyz" << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
        std::ifstream ifs(argv[1]);

        if (!ifs)
            throw std::runtime_error("could not open the XYZ file");

        std::string comment;
        kit::io::load_xyz(ifs, comment, elements, crd);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    if (!kit::is_water(elements, crd, true)) {
        std::cerr << " ** Error ** : not water?" << std::endl;
        return EXIT_FAILURE;
    }

    hline(std::cout);
    const size_t nw = elements.size()/3;
    kit::io::save_xyz(std::cout, argv[1], elements, crd);

#   define DO_ONE(t) \
    { \
        t pot; \
        invoke_potential(pot, nw, &(crd[0])); \
    }

    DO_ONE(h2o::qtip4pf)

    DO_ONE(h2o::ttm2f)
    DO_ONE(h2o::ttm3f)
    DO_ONE(h2o::ttm4f)

    hline(std::cout);

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
