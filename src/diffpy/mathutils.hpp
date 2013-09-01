/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
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

// binary functor for round-off aware comparison

class EpsilonCompare : public std::binary_function<double, double, bool>
{
    public:

        // constructor
        EpsilonCompare(double eps=SQRT_DOUBLE_EPS) : meps(eps)  { }

        // comparison method
        template <class T1, class T2>
        bool operator()(const T1& x, const T2& y)
        {
            return eps_lt(x, y, meps);
        }

    private:

        double meps;
};

// compare if all elements in 2 containers are close.

template<class Seq0, class Seq1>
    bool allclose(const Seq0& v0, const Seq1& v1, double eps=SQRT_DOUBLE_EPS)
{
    bool rv = (v0.size() == v1.size()) &&
        std::equal(v0.begin(), v0.end(), v1.begin(), EpsilonCompare(eps));
    return rv;
}

}   // namespace mathutils
}   // namespace diffpy

// Implementation ------------------------------------------------------------

#include <diffpy/mathutils.ipp>

#endif  // MATHUTILS_HPP_INCLUDED
