/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Implementation of template functions declared in diffpy/serialization.hpp:
*
*   diffpy::serialization_tostring
*   diffpy::serialization_fromstring
*
* Definition of macros for explicit instantiation of serialization_tostring
* and serialization_fromstring for type C or type boost::shared_ptr<C>.
*
*   DIFFPY_INSTANTIATE_SERIALIZATION(C)
*   DIFFPY_INSTANTIATE_PTR_SERIALIZATION(C)
*
*****************************************************************************/

#ifndef SERIALIZATION_IPP_INCLUDED
#define SERIALIZATION_IPP_INCLUDED

#include <sstream>
#include <boost/shared_ptr.hpp>
#include <diffpy/serialization.hpp>

namespace diffpy {

// Template functions --------------------------------------------------------

template <class T>
std::string serialization_tostring(const T& tobj)
{
    std::ostringstream storage(std::ios::binary);
    serialization::oarchive oa(storage, std::ios::binary);
    oa << tobj;
    return storage.str();
}


template <class T>
void serialization_fromstring(T& tobj, const std::string& s)
{
    std::istringstream storage(s, std::ios::binary);
    serialization::iarchive ia(storage, std::ios::binary);
    ia >> tobj;
}

}   // namespace diffpy

// Macros --------------------------------------------------------------------

/// Insert explicit instantiations for serialization_tostring and
/// serialization_fromstring for class C.

#define DIFFPY_INSTANTIATE_SERIALIZATION(C) \
    template \
        std::string \
        diffpy::serialization_tostring<C>(const C&); \
    template \
        void \
        diffpy::serialization_fromstring<C>(C&, const std::string&); \

/// Insert explicit instantiations for serialization_tostring and
/// serialization_fromstring for type boost::shared_ptr<C>.

#define DIFFPY_INSTANTIATE_PTR_SERIALIZATION(C) \
    DIFFPY_INSTANTIATE_SERIALIZATION(boost::shared_ptr<C>)

#endif  // SERIALIZATION_IPP_INCLUDED
