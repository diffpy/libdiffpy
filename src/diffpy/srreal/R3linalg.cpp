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

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include <diffpy/hash.hpp>
#include <diffpy/srreal/R3linalg.hpp>

namespace diffpy {
namespace srreal {
namespace R3 {

// Constants
DLL_EXPORT const int Ndim = 3;

// Functions
const Matrix& identity() {
  static const Matrix mx(1, 0, 0, 0, 1, 0, 0, 0, 1);
  return mx;
}

const Matrix& zeromatrix() {
  static const Matrix mx;
  return mx;
}

namespace {

/// LU decomposition with partial pivoting, P * A = L * U.  Unlike cofactor
/// expansion it stays accurate for nearly singular matrices.
struct LUDecomposition {
  /// U on and above the diagonal, unit lower-triangular L below it.
  double lu[3][3];
  /// Row i of P * A is row perm[i] of A.
  int perm[3] = {0, 1, 2};
  /// Odd number of row exchanges, which flips the sign of determinant.
  bool odd = false;
  /// Some pivot is exactly zero.
  bool singular = false;

  explicit LUDecomposition(const Matrix& A) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) lu[i][j] = A(i, j);
    }
    for (int k = 0; k < 3; ++k) {
      int p = k;
      for (int i = k + 1; i < 3; ++i) {
        if (std::fabs(lu[i][k]) > std::fabs(lu[p][k])) p = i;
      }
      if (p != k) {
        std::swap(lu[p], lu[k]);
        std::swap(perm[p], perm[k]);
        odd = !odd;
      }
      // zero pivot means the rest of the column is zero as well
      if (lu[k][k] == 0.0) {
        singular = true;
        continue;
      }
      for (int i = k + 1; i < 3; ++i) {
        const double f = (lu[i][k] /= lu[k][k]);
        for (int j = k + 1; j < 3; ++j) lu[i][j] -= f * lu[k][j];
      }
    }
  }
};

}  // namespace

double determinant(const Matrix& A) {
  const LUDecomposition d(A);
  const double det = d.lu[0][0] * d.lu[1][1] * d.lu[2][2];
  return d.odd ? -det : det;
}

Matrix inverse(const Matrix& A) {
  const LUDecomposition d(A);
  if (d.singular) throw std::invalid_argument("singular matrix");
  Matrix B;
  // solve L * U * x = P * e_j for every column j of the identity matrix
  for (int j = 0; j < 3; ++j) {
    double x[3];
    for (int i = 0; i < 3; ++i) {
      x[i] = (d.perm[i] == j) ? 1.0 : 0.0;
      for (int k = 0; k < i; ++k) x[i] -= d.lu[i][k] * x[k];
    }
    for (int i = 2; i >= 0; --i) {
      for (int k = i + 1; k < 3; ++k) x[i] -= d.lu[i][k] * x[k];
      x[i] /= d.lu[i][i];
      B(i, j) = x[i];
    }
  }
  return B;
}

size_t hash_value(const Vector& v) {
  return diffpy::hash_range(v.begin(), v.end());
}

size_t hash_value(const Matrix& A) {
  return diffpy::hash_range(A.data().begin(), A.data().end());
}

}  // namespace R3
}  // namespace srreal

namespace mathutils {

// EpsilonLess specialization ------------------------------------------------

template <>
bool EpsilonLess::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
  const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const {
  return std::ranges::lexicographical_compare(A.data(), B.data(), *this);
}

// EpsilonEqual specializations ----------------------------------------------

template <>
bool EpsilonEqual::operator()<srreal::R3::Vector, srreal::R3::Vector>(
  const srreal::R3::Vector& u, const srreal::R3::Vector& v) const {
  return std::ranges::equal(u, v, *this);
}

template <>
bool EpsilonEqual::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
  const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const {
  return std::ranges::equal(A.data(), B.data(), *this);
}

}  // namespace mathutils
}  // namespace diffpy

// End of file
