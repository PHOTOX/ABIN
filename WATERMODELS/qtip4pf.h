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

#ifndef QTIP4PF_H
#define QTIP4PF_H

#include <cstddef>

namespace h2o {

struct qtip4pf {

    qtip4pf();
    ~qtip4pf();

    const char* name() const
    {
        return "q-TIP4P/F";
    }

    double operator()(size_t nw, const double*); // O H H O H H
    double operator()(size_t nw, const double*, double*); // O H H O H H

private:
    size_t  m_nw;
    double* m_mem;

    void allocate(size_t);
};

} // namespace h2o

#endif // QTIP4PF_H
