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
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Declaration of template functions:
*   diffpy::serialization_tostring
*   diffpy::serialization_fromstring
*
* Shared definitions of serialization archive types:
*   diffpy::serialization::iarchive
*   diffpy::serialization::oarchive
*
*****************************************************************************/

#ifndef SERIALIZATION_HPP_INCLUDED
#define SERIALIZATION_HPP_INCLUDED

#include <string>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace diffpy {

template <typename T>
    std::string serialization_tostring(const T& tobj);
template <typename T>
    void serialization_fromstring(T& tobj, const std::string& s);

namespace serialization {

typedef ::boost::archive::binary_iarchive iarchive;
typedef ::boost::archive::binary_oarchive oarchive;

}   // namespace serialization
}   // namespace diffpy

#endif  // SERIALIZATION_HPP_INCLUDED
