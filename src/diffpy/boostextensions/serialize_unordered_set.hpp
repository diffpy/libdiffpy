/*****************************************************************************
*
* Boost serialization support for unordered_set and unordered_multiset.
* Converted from <boost/serialization/hash_set.hpp>
*
* $Id$
*
*****************************************************************************/

#ifndef SERIALIZE_UNORDERED_SET_HPP_INCLUDED
#define SERIALIZE_UNORDERED_SET_HPP_INCLUDED

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// hash_set.hpp: serialization for stl hash_set templates

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>
#include <boost/unordered_set.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/split_free.hpp>

namespace boost { 
namespace serialization {

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

// hash_multiset
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

#include <boost/serialization/collection_traits.hpp>

BOOST_SERIALIZATION_COLLECTION_TRAITS(boost::unordered_set)
BOOST_SERIALIZATION_COLLECTION_TRAITS(boost::unordered_multiset)

#endif  // SERIALIZE_UNORDERED_SET_HPP_INCLUDED
