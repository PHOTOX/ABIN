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

#ifndef IO_XYZ_H
#define IO_XYZ_H

#include <vector>
#include <string>
#include <iostream>

namespace kit { namespace io {

void load_xyz(std::istream&             is,
              std::string&              comment,
              std::vector<std::string>& elements, // clears the original
              std::vector<double>&      xyz);     // clears the original

void save_xyz(std::ostream&                   os,
              const std::string&              comment,
              const std::vector<std::string>& elements,
              const std::vector<double>&      xyz);

}} // namespace kit::io

#endif // IO_XYZ_H
