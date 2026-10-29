/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.rst for a list of people who contributed.
* See LICENSE.rst for license information.
*
******************************************************************************
*
* Implementation of template functions declared in diffpy/serialization.hpp:
*
*   diffpy::serialization_tostring
*   diffpy::serialization_fromstring
*
* Definition of macros for explicit instantiation of serialization_tostring
* and serialization_fromstring for type C or type std::shared_ptr<C>.
*
*   DIFFPY_INSTANTIATE_SERIALIZATION(C)
*   DIFFPY_INSTANTIATE_PTR_SERIALIZATION(C)
*
*****************************************************************************/

#ifndef SERIALIZATION_IPP_INCLUDED
#define SERIALIZATION_IPP_INCLUDED

#include <sstream>
#include <memory>
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

#ifdef _WIN32

// On Windows serialization.hpp includes this file, so the templates are
// instantiated implicitly where used.  Skip the explicit instantiations to
// avoid dllimport/dllexport clashes.
#define DIFFPY_INSTANTIATE_SERIALIZATION(C)
#define DIFFPY_INSTANTIATE_PTR_SERIALIZATION(C)

#else

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
/// serialization_fromstring for type std::shared_ptr<C>.

#define DIFFPY_INSTANTIATE_PTR_SERIALIZATION(C) \
    DIFFPY_INSTANTIATE_SERIALIZATION(std::shared_ptr<C>)

#endif  // _WIN32

#endif  // SERIALIZATION_IPP_INCLUDED
