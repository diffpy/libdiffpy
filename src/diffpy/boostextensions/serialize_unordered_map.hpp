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
* Boost serialization support for unordered_map and unordered_multimap.
* Workaround with conversion to STL list.
*
* $Id$
*
*****************************************************************************/

#ifndef SERIALIZE_UNORDERED_MAP_HPP_INCLUDED
#define SERIALIZE_UNORDERED_MAP_HPP_INCLUDED

#include <boost/config.hpp>
#include <boost/unordered_map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/split_free.hpp>

namespace boost {
namespace serialization {

// unordered_map -------------------------------------------------------------

template<
    class Archive,
    class Key,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void save(
    Archive & ar,
    const boost::unordered_map<
        Key, T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    typedef std::pair<const Key, T> value_type;
    std::list<value_type, Alloc> values(t.begin(), t.end());
    ar & values;
}


template<
    class Archive,
    class Key,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void load(
    Archive & ar,
    boost::unordered_map<
        Key, T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    typedef std::pair<const Key, T> value_type;
    std::list<value_type, Alloc> values;
    ar & values;
    t.clear();
    t.insert(values.begin(), values.end());
}

// split non-intrusive serialization function member into separate
// non intrusive save/load member functions

template<
    class Archive,
    class Key,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void serialize(
    Archive & ar,
    boost::unordered_map<
        Key, T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    boost::serialization::split_free(ar, t, file_version);
}

// unordered_multimap --------------------------------------------------------

template<
    class Archive,
    class Key,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void save(
    Archive & ar,
    const boost::unordered_multimap<
        Key, T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    typedef std::pair<const Key, T> value_type;
    std::list<value_type, Alloc> values(t.begin(), t.end());
    ar & values;
}


template<
    class Archive,
    class Key,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void load(
    Archive & ar,
    boost::unordered_multimap<
        Key, T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    typedef std::pair<const Key, T> value_type;
    std::list<value_type, Alloc> values;
    ar & values;
    t.clear();
    t.insert(values.begin(), values.end());
}

// split non-intrusive serialization function member into separate
// non intrusive save/load member functions

template<
    class Archive,
    class Key,
    class T,
    class Hash,
    class Pred,
    class Alloc
>
inline void serialize(
    Archive & ar,
    boost::unordered_multimap<
        Key, T, Hash, Pred, Alloc
    > &t,
    const unsigned int file_version
){
    boost::serialization::split_free(ar, t, file_version);
}

} // namespace serialization
} // namespace boost

#endif  // SERIALIZE_UNORDERED_MAP_HPP_INCLUDED
