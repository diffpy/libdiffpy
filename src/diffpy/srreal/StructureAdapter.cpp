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
 * class StructureAdapter -- abstract base class for interfacing general
 *     structure objects with srreal classes such as PairQuantity
 *
 *****************************************************************************/

#include <cassert>
#include <cctype>

#include <diffpy/serialization.ipp>
#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/StructureDifference.hpp>

using namespace std;
using diffpy::mathutils::eps_eq;

namespace diffpy {
namespace srreal {

// Public Methods ------------------------------------------------------------

double StructureAdapter::totalOccupancy() const {
  double total_occupancy = 0.0;
  int cntsites = this->countSites();
  for (int i = 0; i < cntsites; ++i) {
    total_occupancy += this->siteOccupancy(i) * this->siteMultiplicity(i);
  }
  return total_occupancy;
}

double StructureAdapter::numberDensity() const { return 0.0; }

const string& StructureAdapter::siteAtomType(int idx) const {
  static string rv = "";
  return rv;
}

int StructureAdapter::siteMultiplicity(int idx) const { return 1; }

double StructureAdapter::siteOccupancy(int idx) const { return 1.0; }

StructureDifference StructureAdapter::diff(
  StructureAdapterConstPtr other) const {
  StructureDifference sd(this->shared_from_this(), other);
  return sd;
}

// Routines ------------------------------------------------------------------

double meanSquareDisplacement(const R3::Matrix& Uijcartn, const R3::Vector& s,
                              bool anisotropy) {
  const R3::Matrix& U = Uijcartn;
  if (!anisotropy) {
    assert(eps_eq(U(0, 0), U(1, 1)));
    assert(eps_eq(U(0, 0), U(2, 2)));
    return U(0, 0);
  }
  assert(R3::norm(s) > 0);
  assert(eps_eq(U(0, 1), U(1, 0)));
  assert(eps_eq(U(0, 2), U(2, 0)));
  assert(eps_eq(U(1, 2), U(2, 1)));
  // s'Us / s's is the MSD along s without normalizing s first.
  const double s0 = s[0], s1 = s[1], s2 = s[2];
  const double sUs =
    U(0, 0) * s0 * s0 + U(1, 1) * s1 * s1 + U(2, 2) * s2 * s2 +
    2 * (U(0, 1) * s0 * s1 + U(0, 2) * s0 * s2 + U(1, 2) * s1 * s2);
  return sUs / (s0 * s0 + s1 * s1 + s2 * s2);
}

double maxUii(StructureAdapterPtr stru) {
  if (!stru) return 0.0;
  double rv = 0.0;
  for (int i = 0; i < stru->countSites(); ++i) {
    const R3::Matrix& U = stru->siteCartesianUij(i);
    for (int k = 0; k < R3::Ndim; k++) {
      if (U(k, k) > rv) rv = U(k, k);
    }
  }
  return rv;
}

}  // namespace srreal
}  // namespace diffpy

DIFFPY_INSTANTIATE_PTR_SERIALIZATION(diffpy::srreal::StructureAdapter)

// End of file
