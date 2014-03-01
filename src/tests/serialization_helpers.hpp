/*****************************************************************************
*
* diffpy.lib        by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* helper template function for testing of serialization.
*
*****************************************************************************/

#ifndef SERIALIZATION_HELPERS_HPP_INCLUDED
#define SERIALIZATION_HELPERS_HPP_INCLUDED

#include <string>
#include <diffpy/serialization.hpp>


/// Return a duplicate object by saving and loading from a string
template <class T>
T dumpandload(const T& src)
{
    using namespace diffpy;
    std::string data = serialization_tostring(src);
    T rv;
    serialization_fromstring(rv, data);
    return rv;
}


#endif  // SERIALIZATION_HELPERS_HPP_INCLUDED
