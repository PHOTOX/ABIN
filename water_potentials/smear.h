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

#ifndef SMEAR_H
#define SMEAR_H

namespace ttm {

struct smear {
    virtual ~smear() {};

    virtual double aCC() const = 0;
    virtual double aCD() const = 0;
    virtual double aDD_inter() const = 0;
    virtual double aDD_intra() const = 0;

    virtual void smear01(const double&, const double&, const double&,
                         double&, double&) const = 0;

    virtual void smear2(const double&, const double&, const double&,
                        double&, double&) const = 0;

    virtual void smear3(const double&, const double&, const double&,
                        double&, double&, double&) const = 0;
};

struct smear_ttm3 : public smear {

    double aCC() const { return 0.175; }
    double aCD() const { return 0.175; }
    double aDD_inter() const { return 0.175; }
    double aDD_intra() const { return 0.175; }

    void smear01(const double&, const double&, const double&,
                 double&, double&) const;

    void smear2(const double&, const double&, const double&,
                double&, double&) const;

    void smear3(const double&, const double&, const double&,
                double&, double&, double&) const;
};

struct smear_ttm2 : public smear_ttm3 {

    double aCC() const { return 0.2; }
    double aCD() const { return 0.2; }
    double aDD_inter() const { return 0.3; }
    double aDD_intra() const { return 0.3; }
};

struct smear_ttm4 : public smear {

    double aCC() const { return 0.4; }
    double aCD() const { return 0.4; }
    double aDD_inter() const { return 0.055; }
    double aDD_intra() const { return 0.626; }

    void smear01(const double&, const double&, const double&,
                 double&, double&) const;

    void smear2(const double&, const double&, const double&,
                double&, double&) const;

    void smear3(const double&, const double&, const double&,
                double&, double&, double&) const;
};

} // namespace ttm

#endif // SMEAR_H
