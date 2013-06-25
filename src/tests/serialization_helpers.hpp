/*****************************************************************************
*
* diffpy.lib        by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Michigan State University.
*                   All rights reserved.
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

#include <sstream>
#include <diffpy/serialization.hpp>


/// Return a duplicate object by saving and loading from a string
template <class T>
T dumpandload(const T& src)
{
    using namespace std;
    stringstream storage(ios::in | ios::out | ios::binary);
    diffpy::serialization::oarchive oa(storage, ios::binary);
    oa << src;
    diffpy::serialization::iarchive ia(storage, ios::binary);
    T rv;
    ia >> rv;
    return rv;
}


#endif  // SERIALIZATION_HELPERS_HPP_INCLUDED
