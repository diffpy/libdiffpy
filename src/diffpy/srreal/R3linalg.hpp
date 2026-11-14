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
 * R3linalg -- vector and matrix types and linear algebra operations in R3 space
 *
 *****************************************************************************/

#ifndef R3LINALG_HPP_INCLUDED
#define R3LINALG_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <boost/serialization/access.hpp>
#include <diffpy/mathutils.hpp>

#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {
namespace R3 {

// Constants -----------------------------------------------------------------

extern DLL_EXPORT const int Ndim;
using ::diffpy::mathutils::SQRT_DOUBLE_EPS;

// Classes -------------------------------------------------------------------

class Vector {
 public:
  typedef double value_type;
  typedef std::size_t size_type;
  typedef std::array<double, 3> array_type;
  typedef array_type::iterator iterator;
  typedef array_type::const_iterator const_iterator;

  // constructors
  constexpr Vector() noexcept : mdata{} {}
  constexpr Vector(double x, double y, double z) noexcept : mdata{x, y, z} {}
  constexpr Vector(const Vector&) noexcept = default;
  constexpr Vector(Vector&&) noexcept = default;

  // Assignments are restricted to lvalues so that writes to a temporary,
  // such as row(A, 0) = v, fail to compile instead of being lost.
  constexpr Vector& operator=(const Vector&) & noexcept = default;
  constexpr Vector& operator=(Vector&&) & noexcept = default;

  // Element writes need an lvalue too, row(A, 0)[1] = x does not compile.
  constexpr double& operator[](size_type i) & noexcept { return mdata[i]; }
  constexpr const double& operator[](size_type i) const& noexcept {
    return mdata[i];
  }
  constexpr double& operator()(size_type i) & noexcept { return mdata[i]; }
  constexpr const double& operator()(size_type i) const& noexcept {
    return mdata[i];
  }
  static constexpr size_type size() noexcept { return 3; }
  constexpr array_type& data() & noexcept { return mdata; }
  constexpr const array_type& data() const& noexcept { return mdata; }

  // iterators, writable only for lvalues
  constexpr iterator begin() & noexcept { return mdata.begin(); }
  constexpr iterator end() & noexcept { return mdata.end(); }
  constexpr const_iterator begin() const& noexcept { return mdata.begin(); }
  constexpr const_iterator end() const& noexcept { return mdata.end(); }

  // exact element-wise equality, != is rewritten from it
  constexpr bool operator==(const Vector&) const noexcept = default;

  // in-place arithmetic
  constexpr Vector& operator+=(const Vector& v) & noexcept {
    for (size_type k = 0; k < 3; ++k) mdata[k] += v.mdata[k];
    return *this;
  }

  constexpr Vector& operator-=(const Vector& v) & noexcept {
    for (size_type k = 0; k < 3; ++k) mdata[k] -= v.mdata[k];
    return *this;
  }

  constexpr Vector& operator*=(double s) & noexcept {
    for (double& x : mdata) x *= s;
    return *this;
  }

  constexpr Vector& operator/=(double s) & noexcept {
    for (double& x : mdata) x /= s;
    return *this;
  }

 private:
  array_type mdata;

  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);
};

class Matrix {
 public:
  typedef double value_type;
  typedef std::size_t size_type;
  typedef std::array<double, 9> array_type;

  // constructors
  constexpr Matrix() noexcept : mdata{} {}
  constexpr Matrix(double x0, double x1, double x2, double x3, double x4,
                   double x5, double x6, double x7, double x8) noexcept
      : mdata{x0, x1, x2, x3, x4, x5, x6, x7, x8} {}
  constexpr Matrix(const Matrix&) noexcept = default;
  constexpr Matrix(Matrix&&) noexcept = default;

  // assignments are restricted to lvalues, see Vector
  constexpr Matrix& operator=(const Matrix&) & noexcept = default;
  constexpr Matrix& operator=(Matrix&&) & noexcept = default;

  // row-major element access, writable only for lvalues
  constexpr double& operator()(size_type i, size_type j) & noexcept {
    return mdata[3 * i + j];
  }
  constexpr const double& operator()(size_type i, size_type j) const& noexcept {
    return mdata[3 * i + j];
  }
  static constexpr size_type size1() noexcept { return 3; }
  static constexpr size_type size2() noexcept { return 3; }
  constexpr array_type& data() & noexcept { return mdata; }
  constexpr const array_type& data() const& noexcept { return mdata; }

  // exact element-wise equality, != is rewritten from it
  constexpr bool operator==(const Matrix&) const noexcept = default;

  // in-place arithmetic
  constexpr Matrix& operator+=(const Matrix& A) & noexcept {
    for (size_type k = 0; k < 9; ++k) mdata[k] += A.mdata[k];
    return *this;
  }

  constexpr Matrix& operator-=(const Matrix& A) & noexcept {
    for (size_type k = 0; k < 9; ++k) mdata[k] -= A.mdata[k];
    return *this;
  }

  constexpr Matrix& operator*=(double s) & noexcept {
    for (double& x : mdata) x *= s;
    return *this;
  }

  constexpr Matrix& operator/=(double s) & noexcept {
    for (double& x : mdata) x /= s;
    return *this;
  }

 private:
  array_type mdata;

  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);
};

inline constexpr Vector zerovector{};

// Functions -----------------------------------------------------------------

DLL_EXPORT const Matrix& identity();
DLL_EXPORT const Matrix& zeromatrix();
DLL_EXPORT double determinant(const Matrix& A);
DLL_EXPORT Matrix inverse(const Matrix& A);

// Vector arithmetic ---------------------------------------------------------

constexpr Vector operator-(const Vector& u) noexcept {
  return Vector(-u[0], -u[1], -u[2]);
}

constexpr Vector operator+(Vector u, const Vector& v) noexcept {
  return u += v;
}

constexpr Vector operator-(Vector u, const Vector& v) noexcept {
  return u -= v;
}

constexpr Vector operator*(Vector u, double s) noexcept { return u *= s; }
constexpr Vector operator*(double s, Vector u) noexcept { return u *= s; }
constexpr Vector operator/(Vector u, double s) noexcept { return u /= s; }

// Matrix arithmetic ---------------------------------------------------------

constexpr Matrix operator-(Matrix A) noexcept { return A *= -1.0; }

constexpr Matrix operator+(Matrix A, const Matrix& B) noexcept {
  return A += B;
}

constexpr Matrix operator-(Matrix A, const Matrix& B) noexcept {
  return A -= B;
}

constexpr Matrix operator*(Matrix A, double s) noexcept { return A *= s; }
constexpr Matrix operator*(double s, Matrix A) noexcept { return A *= s; }
constexpr Matrix operator/(Matrix A, double s) noexcept { return A /= s; }

/// Matrix product A * B.
constexpr Matrix prod(const Matrix& A, const Matrix& B) noexcept {
  Matrix C;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      C(i, j) = A(i, 0) * B(0, j) + A(i, 1) * B(1, j) + A(i, 2) * B(2, j);
    }
  }
  return C;
}

constexpr Matrix trans(const Matrix& A) noexcept {
  return Matrix(A(0, 0), A(1, 0), A(2, 0), A(0, 1), A(1, 1), A(2, 1), A(0, 2),
                A(1, 2), A(2, 2));
}

/// Copy of the i-th row of A.
constexpr Vector row(const Matrix& A, std::size_t i) noexcept {
  return Vector(A(i, 0), A(i, 1), A(i, 2));
}

/// Copy of the j-th column of A.
constexpr Vector column(const Matrix& A, std::size_t j) noexcept {
  return Vector(A(0, j), A(1, j), A(2, j));
}

constexpr void setrow(Matrix& A, std::size_t i, const Vector& v) noexcept {
  A(i, 0) = v[0];
  A(i, 1) = v[1];
  A(i, 2) = v[2];
}

constexpr void setcolumn(Matrix& A, std::size_t j, const Vector& v) noexcept {
  A(0, j) = v[0];
  A(1, j) = v[1];
  A(2, j) = v[2];
}

// Hashing -------------------------------------------------------------------

DLL_EXPORT size_t hash_value(const Vector& v);
DLL_EXPORT size_t hash_value(const Matrix& A);

// Vector functions ----------------------------------------------------------

inline Vector floor(const Vector& v) {
  return Vector(std::floor(v[0]), std::floor(v[1]), std::floor(v[2]));
}

template <class V>
double dot(const V& u, const V& v) {
  return (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}

template <class V>
double norm(const V& u) {
  return std::sqrt(R3::dot(u, u));
}

template <class V>
double distance(const V& u, const V& v) {
  const double dx = u[0] - v[0];
  const double dy = u[1] - v[1];
  const double dz = u[2] - v[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

template <class V>
Vector cross(const V& u, const V& v) {
  return Vector(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0]);
}

/// Matrix-vector product M * u.
template <class V>
Vector mxvecproduct(const Matrix& M, const V& u) {
  return Vector(M(0, 0) * u[0] + M(0, 1) * u[1] + M(0, 2) * u[2],
                M(1, 0) * u[0] + M(1, 1) * u[1] + M(1, 2) * u[2],
                M(2, 0) * u[0] + M(2, 1) * u[1] + M(2, 2) * u[2]);
}

/// Vector-matrix product u * M.
template <class V>
Vector mxvecproduct(const V& u, const Matrix& M) {
  return Vector(u[0] * M(0, 0) + u[1] * M(1, 0) + u[2] * M(2, 0),
                u[0] * M(0, 1) + u[1] * M(1, 1) + u[2] * M(2, 1),
                u[0] * M(0, 2) + u[1] * M(1, 2) + u[2] * M(2, 2));
}

}  // namespace R3
}  // namespace srreal

namespace mathutils {

// Template specializations --------------------------------------------------

template <>
DLL_EXPORT bool EpsilonLess::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
  const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const;

template <>
DLL_EXPORT bool
EpsilonEqual::operator()<srreal::R3::Vector, srreal::R3::Vector>(
  const srreal::R3::Vector& u, const srreal::R3::Vector& v) const;

template <>
DLL_EXPORT bool
EpsilonEqual::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
  const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const;

}  // namespace mathutils
}  // namespace diffpy

// Implementation of serialization
#include <diffpy/srreal/R3linalg.ipp>

#endif  // R3LINALG_HPP_INCLUDED
