/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * Hash utilities that replace Boost.Functional/Hash:
 *
 *   diffpy::hash_of      -- hash_value(v) found by ADL, else std::hash<T>
 *   diffpy::hash_combine -- mix the hash of a value into a seed
 *   diffpy::hash_range   -- hash of a sequence of values
 *   diffpy::hash_values  -- hash of several values
 *   diffpy::PairHash     -- hash function object for std::pair
 *
 *****************************************************************************/

#ifndef HASH_HPP_INCLUDED
#define HASH_HPP_INCLUDED

#include <concepts>
#include <cstddef>
#include <functional>
#include <utility>

namespace diffpy {

namespace detail {

/// True when hash_value(v) can be found by argument-dependent lookup.
template <class T>
concept HasHashValue = requires(const T& v) {
  { hash_value(v) } -> std::convertible_to<std::size_t>;
};

}  // namespace detail

/// Hash of v from hash_value(v) when available, otherwise std::hash<T>.
template <class T>
std::size_t hash_of(const T& v) {
  if constexpr (detail::HasHashValue<T>) {
    return hash_value(v);
  } else {
    return std::hash<T>{}(v);
  }
}

/// Mix the hash of v into seed, same scheme as boost::hash_combine.
template <class T>
void hash_combine(std::size_t& seed, const T& v) {
  const std::size_t magic = static_cast<std::size_t>(0x9e3779b97f4a7c15ull);
  seed ^= hash_of(v) + magic + (seed << 6) + (seed >> 2);
}

/// Order-dependent hash of the values in [first, last).
template <class It>
std::size_t hash_range(It first, It last) {
  std::size_t seed = 0;
  for (; first != last; ++first) hash_combine(seed, *first);
  return seed;
}

/// Order-dependent hash of the arguments, same as hash_combine of each.
template <class... Ts>
std::size_t hash_values(const Ts&... values) {
  std::size_t seed = 0;
  (hash_combine(seed, values), ...);
  return seed;
}

/// Hash function object for std::pair, for use in unordered containers.
struct PairHash {
  template <class A, class B>
  std::size_t operator()(const std::pair<A, B>& p) const {
    return hash_values(p.first, p.second);
  }
};

}  // namespace diffpy

#endif  // HASH_HPP_INCLUDED
