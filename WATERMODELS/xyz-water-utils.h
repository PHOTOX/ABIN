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

#ifndef XYZ_WATER_UTILS_H
#define XYZ_WATER_UTILS_H

#include <string>
#include <vector>

namespace kit {

bool is_water(std::vector<std::string>& elements,
              std::vector<double>&      xyz,
              bool reorder = true, const double& rOH = 1.3);

void make_water_elements(size_t nw, std::vector<std::string>& elements);

} // namespace kit

#endif // XYZ_WATER_UTILS_H
