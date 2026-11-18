/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 *****************************************************************************/

#ifndef PAIR_TEST_HELPERS_HPP_INCLUDED
#define PAIR_TEST_HELPERS_HPP_INCLUDED

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <diffpy/srreal/NoMetaStructureAdapter.hpp>
#include <diffpy/srreal/NoSymmetryStructureAdapter.hpp>
#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/srreal/ScatteringFactorTable.hpp>

#include "test_helpers.hpp"

namespace pairtests {
using namespace diffpy::srreal;

using Setup = std::function<void(PairQuantity&)>;

inline std::shared_ptr<PeriodicStructureAdapter> loadPeriodic(
  const std::string& file, bool varyOccupancy = false) {
  auto rv = std::dynamic_pointer_cast<PeriodicStructureAdapter>(
    loadTestPeriodicStructure(file));
  REQUIRE(rv);
  for (int i = 0; i < rv->countSites(); ++i) {
    Atom& a = rv->at(i);
    if (a.uij_cartn == R3::zeromatrix()) {
      a.uij_cartn = (0.004 + 0.001 * (i % 3)) * R3::identity();
    }
    if (varyOccupancy) a.occupancy = 1.0 - 0.1 * (i % 3);
  }
  return rv;
}

inline std::shared_ptr<AtomicStructureAdapter> asCluster(
  const PeriodicStructureAdapter& stru) {
  auto rv = std::make_shared<AtomicStructureAdapter>();
  rv->assign(stru.begin(), stru.end());
  return rv;
}

/// The first half of the atoms with the operations of a space group that
/// suits the lattice.
inline std::shared_ptr<CrystalStructureAdapter> asCrystal(
  const PeriodicStructureAdapter& stru) {
  const Lattice& L = stru.getLattice();
  auto rv = std::make_shared<CrystalStructureAdapter>();
  rv->setLatPar(L.a(), L.b(), L.c(), L.alpha(), L.beta(), L.gamma());
  const int n = std::max(1, stru.countSites() / 2);
  for (int i = 0; i < n; ++i) rv->append(stru[i]);
  rv->addSymOp(R3::identity(), R3::zerovector);
  rv->addSymOp(-1.0 * R3::identity(), R3::zerovector);
  if (L.alpha() == 90.0 && L.beta() == 90.0) {
    R3::Matrix R = R3::identity();
    R(0, 0) = -1.0;
    R(1, 1) = -1.0;
    rv->addSymOp(R, R3::Vector(0.5, 0.5, 0.0));
    rv->addSymOp(-1.0 * R, R3::Vector(0.5, 0.5, 0.0));
  }
  return rv;
}

inline std::vector<std::pair<std::string, StructureAdapterPtr>> allStructures(
  bool varyOccupancy = false) {
  std::vector<std::pair<std::string, StructureAdapterPtr>> rv;
  for (const char* f : {"ZnS_wurtzite.stru", "CaTiO3.stru", "Ni.stru",
                        "alpha_K2Bi8Se13.stru"}) {
    const auto periodic = loadPeriodic(f, varyOccupancy);
    const auto crystal = asCrystal(*periodic);
    const std::string name = f;
    rv.emplace_back(name + " periodic", periodic);
    rv.emplace_back(name + " crystal", crystal);
    rv.emplace_back(name + " cluster", asCluster(*periodic));
    rv.emplace_back(name + " nometa periodic", nometa(periodic));
    rv.emplace_back(name + " nosymmetry crystal", nosymmetry(crystal));
    rv.emplace_back(name + " nometa nosymmetry periodic",
                    nometa(nosymmetry(periodic)));
  }
  return rv;
}

inline void setupDefault(PairQuantity& calc) {
  calc.setRmax(9.5);
  calc.setDoubleAttr("qmax", 22.0);
}

template <class Check>
void forEachCalculatorCase(Check check) {
  for (const auto& [name, structure] : allStructures(true)) {
    INFO(name);
    check(setupDefault, structure);
  }
  const auto perovskite = loadPeriodic("CaTiO3.stru", true);
  const std::pair<const char*, Setup> setups[] = {
    {"jeong",
     [](PairQuantity& c) {
       setupDefault(c);
       c.setDoubleAttr("delta1", 0.3);
       c.setDoubleAttr("delta2", 1.2);
       c.setDoubleAttr("qbroad", 0.03);
       c.setDoubleAttr("qbroad_seperable", 0.02);
     }},
    {"debye-waller",
     [](PairQuantity& c) {
       setupDefault(c);
       dynamic_cast<PeakWidthModelOwner&>(c).setPeakWidthModelByType(
         "debye-waller");
     }},
    {"constant",
     [](PairQuantity& c) {
       setupDefault(c);
       dynamic_cast<PeakWidthModelOwner&>(c).setPeakWidthModelByType(
         "constant");
       c.setDoubleAttr("width", 0.21);
     }},
    {"sharpened peaks",
     [](PairQuantity& c) {
       setupDefault(c);
       c.setDoubleAttr("delta2", 5.5);
     }},
    {"rmin",
     [](PairQuantity& c) {
       setupDefault(c);
       c.setRmin(3.0);
       c.setRmax(11.0);
     }},
    {"neutrons",
     [](PairQuantity& c) {
       setupDefault(c);
       dynamic_cast<ScatteringFactorTableOwner&>(c)
         .setScatteringFactorTableByType("neutron");
     }},
  };
  for (const auto& [name, setup] : setups) {
    INFO(name);
    check(setup, perovskite);
  }
}

}  // namespace pairtests
#endif
