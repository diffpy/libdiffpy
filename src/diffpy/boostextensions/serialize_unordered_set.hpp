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
* Boost serialization support for unordered_set and unordered_multiset.
* Workaround with conversion to STL list.
*
* $Id$
*
*****************************************************************************/

#ifndef SERIALIZE_UNORDERED_SET_HPP_INCLUDED
#define SERIALIZE_UNORDERED_SET_HPP_INCLUDED

#include <boost/config.hpp>
#include <boost/unordered_set.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/split_free.hpp>

namespace boost {
namespace serialization {

// unordered_set -------------------------------------------------------------

template<
    class Archive,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void save(
    Archive & ar,
    const boost::unordered_set<
        T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    std::list<T, Alloc> values(t.begin(), t.end());
    ar & values;
}


template<
    class Archive,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void load(
    Archive & ar,
    boost::unordered_set<
        T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    std::list<T, Alloc> values;
    ar & values;
    t.clear();
    t.insert(values.begin(), values.end());
}

// split non-intrusive serialization function member into separate
// non intrusive save/load member functions

template<
    class Archive,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void serialize(
    Archive & ar,
    boost::unordered_set<
        T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    boost::serialization::split_free(ar, t, file_version);
}

// unordered_multiset --------------------------------------------------------

template<
    class Archive,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void save(
    Archive & ar,
    const boost::unordered_multiset<
        T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    std::list<T, Alloc> values(t.begin(), t.end());
    ar & values;
}


template<
    class Archive,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void load(
    Archive & ar,
    boost::unordered_multiset<
        T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    std::list<T, Alloc> values;
    ar & values;
    t.clear();
    t.insert(values.begin(), values.end());
}

// split non-intrusive serialization function member into separate
// non intrusive save/load member functions

template<
    class Archive,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void serialize(
    Archive & ar,
    boost::unordered_multiset<
        T, Hash, Pred, Alloc
    > & t,
    const unsigned int file_version
){
    boost::serialization::split_free(ar, t, file_version);
}

} // namespace serialization
} // namespace boost

#endif  // SERIALIZE_UNORDERED_SET_HPP_INCLUDED
