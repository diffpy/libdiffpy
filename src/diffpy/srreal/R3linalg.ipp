/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2009 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * R3linalg -- serialization of the Vector and Matrix types.
 *
 *****************************************************************************/

#ifndef R3LINALG_IPP_INCLUDED
#define R3LINALG_IPP_INCLUDED

#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/collection_size_type.hpp>
#include <boost/serialization/nvp.hpp>

namespace diffpy {
namespace srreal {
namespace R3 {
namespace detail {

// These wrappers keep the archive layout of the former ublas base types, so
// older archives and pickles still load.  Each stands for one ublas type and
// must stay a distinct class with default serialization traits, because the
// archive writes class information once per type.

inline void check_legacy_size(bool valid) {
  if (!valid) {
    throw boost::archive::archive_exception(
      boost::archive::archive_exception::array_size_too_short);
  }
}

/// Stands for ublas::bounded_array<double, N>.
template <std::size_t N>
struct LegacyBoundedArray {
  std::array<double, N>& a;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int) {
    boost::serialization::collection_size_type s(N);
    ar& boost::serialization::make_nvp("size", s);
    check_legacy_size(s == N);
    ar& boost::serialization::make_array(a.data(), N);
  }
};

/// Stands for ublas::vector<double, ublas::bounded_array<double, 3> >.
struct LegacyVectorBase {
  std::array<double, 3>& a;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int) {
    LegacyBoundedArray<3> data{a};
    ar& boost::serialization::make_nvp("data", data);
  }
};

/// Stands for ublas::matrix<double, row_major, bounded_array<double, 9> >.
struct LegacyMatrixBase {
  std::array<double, 9>& a;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int) {
    boost::serialization::collection_size_type s1(3), s2(3);
    ar& boost::serialization::make_nvp("size1", s1);
    ar& boost::serialization::make_nvp("size2", s2);
    check_legacy_size(s1 == 3 && s2 == 3);
    LegacyBoundedArray<9> data{a};
    ar& boost::serialization::make_nvp("data", data);
  }
};

}  // namespace detail

template <class Archive>
void Vector::serialize(Archive& ar, const unsigned int) {
  detail::LegacyVectorBase base{mdata};
  ar & base;
}

template <class Archive>
void Matrix::serialize(Archive& ar, const unsigned int) {
  detail::LegacyMatrixBase base{mdata};
  ar & base;
}

}  // namespace R3
}  // namespace srreal
}  // namespace diffpy

#endif  // R3LINALG_IPP_INCLUDED
