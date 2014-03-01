/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* Convenience functions for argument checking.
*
*****************************************************************************/

#ifndef VALIDATORS_HPP_INCLUDED
#define VALIDATORS_HPP_INCLUDED

#include <stdexcept>
#include <string>

#include <diffpy/mathutils.hpp>

namespace diffpy {
namespace validators {

/// Throw invalid_argument for a negative value

template <class T>
void ensureNonNegative(const std::string& vname, const T& value)
{
    if (value < 0)
    {
        std::string emsg(vname);
        emsg += " cannot be negative.";
        throw std::invalid_argument(emsg);
    }
}

/// Throw invalid_argument for a value that is not greater than DOUBLE_EPS

template <class T>
void ensureEpsilonPositive(const std::string& vname, const T& value)
{
    using diffpy::mathutils::eps_gt;
    if (!eps_gt(double(value), 0.0))
    {
        std::string emsg(vname);
        emsg += " must be epsilon positive.";
        throw std::invalid_argument(emsg);
    }
}

/// Throw invalid_argument if the argument is not true.  For boost smart
/// pointers this is equivalent to a truth check of p.get().

template <class T>
void ensureNonNull(const std::string& vname, const T& p)
{
    if (!p)
    {
        std::string emsg(vname);
        emsg += " cannot be NULL.";
        throw std::invalid_argument(emsg);
    }
}

/// Throw runtime_error if a file argument is not true.

template <class T>
void ensureFileOK(const std::string& fname, const T& fp)
{
    if (!fp)
    {
        std::string emsg = "Cannot open '";
        emsg += fname + "'.";
        throw std::runtime_error(emsg);
    }
}

}   // namespace validators
}   // namespace diffpy

#endif  // VALIDATORS_HPP_INCLUDED
