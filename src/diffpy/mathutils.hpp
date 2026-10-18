/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* Various common mathematical constants and functions.
*
*****************************************************************************/

#ifndef MATHUTILS_HPP_INCLUDED
#define MATHUTILS_HPP_INCLUDED

#include <limits>
#include <cmath>
#include <functional>
#include <algorithm>

namespace diffpy {
namespace mathutils {

// constants

const double DOUBLE_MAX = std::numeric_limits<double>().max();
const double DOUBLE_EPS = std::numeric_limits<double>().epsilon();
const double SQRT_DOUBLE_EPS = (sqrt(DOUBLE_EPS) + 1.0) - 1.0;
const double GAUSS_SIGMA_TO_FWHM = 2 * sqrt(2 * M_LN2);

// trigonometric functions with more exact values at n*30 degrees

double cosd(double);
double sind(double x);
double acosd(double x);
double asind(double x);

// round-off aware comparison operations

bool eps_eq(const double& x, const double& y, double eps=SQRT_DOUBLE_EPS);
bool eps_gt(const double& x, const double& y, double eps=SQRT_DOUBLE_EPS);
bool eps_lt(const double& x, const double& y, double eps=SQRT_DOUBLE_EPS);

// binary functor for round-off aware less-than comparison

class EpsilonLess
{
    public:

        // constructor
        EpsilonLess(double eps=SQRT_DOUBLE_EPS) : meps(eps)  { }

        // comparison method
        bool operator()(const double& x, const double& y) const
        {
            return eps_lt(x, y, meps);
        }
        template<class Seq0, class Seq1>
        bool operator()(const Seq0& v0, const Seq1& v1) const
        {
            bool rv = std::lexicographical_compare(
                    v0.begin(), v0.end(),
                    v1.begin(), v1.end(), *this);
            return rv;
        }

    private:

        double meps;
};

// binary functor for round-off aware equality comparison

class EpsilonEqual
{
    public:

        // constructor
        EpsilonEqual(double eps=SQRT_DOUBLE_EPS) : meps(eps)  { }

        // comparison method
        bool operator()(const double& x, const double& y) const
        {
            return eps_eq(x, y, meps);
        }
        template<class Seq0, class Seq1>
        bool operator()(const Seq0& v0, const Seq1& v1) const
        {
            bool rv = (v0.size() == v1.size()) &&
                std::equal(v0.begin(), v0.end(), v1.begin(), *this);
                return rv;
        }

    private:

        double meps;
};


}   // namespace mathutils
}   // namespace diffpy

// Implementation ------------------------------------------------------------

#include <diffpy/mathutils.ipp>

#endif  // MATHUTILS_HPP_INCLUDED
