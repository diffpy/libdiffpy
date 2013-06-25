/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
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
* Shared definitions of serialization archive types:
*   diffpy::serialization::iarchive
*   diffpy::serialization::oarchive
*
*****************************************************************************/

#ifndef SERIALIZATION_HPP_INCLUDED
#define SERIALIZATION_HPP_INCLUDED

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

/// this macro provides to explicit template instantiation of the serialize
/// method for a specified class.
#define DIFFPY_INSTANTIATE_SERIALIZE(C) \
    template void C::serialize( \
            ::diffpy::serialization::iarchive&, const unsigned int); \
    template void C::serialize( \
            ::diffpy::serialization::oarchive&, const unsigned int); \

namespace diffpy {
namespace serialization {

typedef ::boost::archive::binary_iarchive iarchive;
typedef ::boost::archive::binary_oarchive oarchive;

}   // namespace serialization
}   // namespace diffpy

#endif  // SERIALIZATION_HPP_INCLUDED
