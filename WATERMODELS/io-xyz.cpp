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

#include <cassert>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "io-xyz.h"

namespace kit { namespace io {

void save_xyz(std::ostream&                   os,
              const std::string&              comment,
              const std::vector<std::string>& elements,
              const std::vector<double>&      xyz)
{
    assert(3*elements.size() == xyz.size());

    std::ios::fmtflags saved_flags = os.flags();

    os << std::setw(5) << elements.size()
       << '\n' << comment << '\n'
       << std::setprecision(9) << std::scientific;

    for (size_t n = 0; n < elements.size(); ++n)
        os << std::setw(5) << std::left << elements[n]
           << std::setw(18) << std::right << xyz[3*n + 0]
           << std::setw(18) << std::right << xyz[3*n + 1]
           << std::setw(18) << std::right << xyz[3*n + 2]
           << '\n';

    os.flags(saved_flags);
}

void load_xyz(std::istream&             is,
              std::string&              comment,
              std::vector<std::string>& elements,
              std::vector<double>&      xyz)
{
    comment.clear();
    elements.clear();
    xyz.clear();

    std::string line;
    std::getline(is, line);

    if (is.eof())
        return;

    size_t natoms;
    std::istringstream iss(line);
    iss >> natoms;

    xyz.reserve(3*natoms);
    elements.reserve(natoms);

    std::getline(is, comment);

    size_t lineno(2), n(0);
    while (!is.eof()) {
        std::getline(is, line);
        ++lineno;

        if (line.length() == 0)
            continue;

        std::string element;
        double x, y, z;

        iss.clear();
        iss.str(line);
        iss >> element >> x >> y >> z;
        if (iss.fail()) {
            std::ostringstream oss;
            oss << "unexpected text at line " << lineno << " of the XYZ stream";
            throw std::runtime_error(oss.str());
        }

        elements.push_back(element);

        xyz.push_back(x);
        xyz.push_back(y);
        xyz.push_back(z);

        if (++n == natoms)
            break;
    }

    if (n != natoms) {
        std::ostringstream oss;
        oss << "wrong number of atoms ("
            << n << " instead of " << natoms
            << ") in the XYZ stream";
        throw std::runtime_error(oss.str());
    }
}

}} // namespace kit::io
