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
 * TestPeriodicStructureAdapter -- unit tests for an adapter
 *     for periodic structures
 *
 * TestPeriodicStructureBondGenerator -- unit tests for bond generator
 *
 *****************************************************************************/

#include <typeinfo>
#include <sstream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include "test_helpers.hpp"
#include "serialization_helpers.hpp"

using Catch::Matchers::WithinAbs;

namespace diffpy {
namespace srreal {

using namespace std;

// Local Helpers -------------------------------------------------------------

namespace {

int countBonds(BaseBondGenerator& bnds) {
  int rv = 0;
  for (bnds.rewind(); !bnds.finished(); bnds.next()) ++rv;
  return rv;
}

template <class Tstru, class Tbnds>
double testmsd0(const Tstru& stru, const Tbnds& bnds) {
  bool anisotropy0 = stru->siteAnisotropy(bnds->site0());
  double rv =
    meanSquareDisplacement(bnds->Ucartesian0(), bnds->r01(), anisotropy0);
  return rv;
}

template <class Tstru, class Tbnds>
double testmsd1(const Tstru& stru, const Tbnds& bnds) {
  bool anisotropy1 = stru->siteAnisotropy(bnds->site1());
  double rv =
    meanSquareDisplacement(bnds->Ucartesian1(), bnds->r01(), anisotropy1);
  return rv;
}

}  // namespace

//////////////////////////////////////////////////////////////////////////////
// TestPeriodicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestPeriodicStructureAdapterFixture {
 protected:
  StructureAdapterPtr m_ni;
  StructureAdapterPtr m_kbise;
  StructureAdapterPtr m_catio3;
  StructureAdapterPtr m_pswt;

 public:
  TestPeriodicStructureAdapterFixture() {
    m_ni = loadTestPeriodicStructure("Ni.stru");
    m_kbise = loadTestPeriodicStructure("alpha_K2Bi8Se13.stru");
    m_catio3 = loadTestPeriodicStructure("CaTiO3.stru");
    m_pswt = loadTestPeriodicStructure("PbScW25TiO3.stru");
  }
};  // class TestPeriodicStructureAdapterFixture

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_typeid",
                 "[PeriodicStructureAdapter]") {
  StructureAdapter& r_ni = *m_ni;
  CHECK(typeid(PeriodicStructureAdapter) == typeid(r_ni));
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_countSites",
                 "[PeriodicStructureAdapter]") {
  CHECK(4 == m_ni->countSites());
  CHECK(23 == m_kbise->countSites());
  CHECK(20 == m_catio3->countSites());
  CHECK(56 == m_pswt->countSites());
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_totalOccupancy",
                 "[PeriodicStructureAdapter]") {
  const double eps = 10 * diffpy::mathutils::DOUBLE_EPS;
  CHECK(4.0 == m_ni->totalOccupancy());
  CHECK(23.0 == m_kbise->totalOccupancy());
  CHECK(20.0 == m_catio3->totalOccupancy());
  CHECK_THAT(40.0, WithinAbs(m_pswt->totalOccupancy(), eps));
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_numberDensity",
                 "[PeriodicStructureAdapter]") {
  const double eps = 1.0e-7;
  CHECK_THAT(0.0914114, WithinAbs(m_ni->numberDensity(), eps));
  CHECK_THAT(0.0335565, WithinAbs(m_kbise->numberDensity(), eps));
  CHECK_THAT(0.0894566, WithinAbs(m_catio3->numberDensity(), eps));
  CHECK_THAT(0.0760332, WithinAbs(m_pswt->numberDensity(), eps));
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_siteCartesianPosition",
                 "[PeriodicStructureAdapter]") {
  const double eps = 1.0e-5;
  R3::Vector rCa = m_catio3->siteCartesianPosition(4);
  CHECK_THAT(0.03486, WithinAbs(rCa[0], eps));
  CHECK_THAT(0.19366, WithinAbs(rCa[1], eps));
  CHECK_THAT(1.90975, WithinAbs(rCa[2], eps));
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_siteAnisotropy",
                 "[PeriodicStructureAdapter]") {
  for (int i = 0; i < m_ni->countSites(); ++i) {
    CHECK(false == m_ni->siteAnisotropy(i));
  }
  for (int i = 0; i < m_catio3->countSites(); ++i) {
    CHECK(true == m_catio3->siteAnisotropy(i));
  }
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_siteCartesianUij",
                 "[PeriodicStructureAdapter]") {
  // nickel should have all Uij equal zero.
  CHECK(R3::zeromatrix() == m_ni->siteCartesianUij(0));
  // check CaTiO3 values
  const R3::Matrix& UO2 = m_catio3->siteCartesianUij(12);
  const double eps = 1e-8;
  CHECK_THAT(0.0065, WithinAbs(UO2(0, 0), eps));
  CHECK_THAT(0.0060, WithinAbs(UO2(1, 1), eps));
  CHECK_THAT(0.0095, WithinAbs(UO2(2, 2), eps));
  CHECK_THAT(0.0020, WithinAbs(UO2(0, 1), eps));
  CHECK_THAT(0.0020, WithinAbs(UO2(1, 0), eps));
  CHECK_THAT(-0.0008, WithinAbs(UO2(0, 2), eps));
  CHECK_THAT(-0.0008, WithinAbs(UO2(2, 0), eps));
  CHECK_THAT(-0.0010, WithinAbs(UO2(1, 2), eps));
  CHECK_THAT(-0.0010, WithinAbs(UO2(2, 1), eps));
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_siteAtomType",
                 "[PeriodicStructureAdapter]") {
  CHECK(string("Ni") == m_ni->siteAtomType(0));
  CHECK(string("Ni") == m_ni->siteAtomType(3));
  CHECK(string("K1+") == m_kbise->siteAtomType(0));
  CHECK(string("Bi3+") == m_kbise->siteAtomType(2));
  CHECK(string("Se") == m_kbise->siteAtomType(10));
  CHECK(string("Se") == m_kbise->siteAtomType(22));
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_getLattice",
                 "[PeriodicStructureAdapter]") {
  PeriodicStructureAdapterPtr pkbise =
    std::dynamic_pointer_cast<PeriodicStructureAdapter>(m_kbise);
  CHECK(pkbise);
  const Lattice& L = pkbise->getLattice();
  const double eps = 1.0e-12;
  CHECK_THAT(13.768, WithinAbs(L.a(), eps));
  CHECK_THAT(12.096, WithinAbs(L.b(), eps));
  CHECK_THAT(4.1656, WithinAbs(L.c(), eps));
  CHECK_THAT(89.98, WithinAbs(L.alpha(), eps));
  CHECK_THAT(98.64, WithinAbs(L.beta(), eps));
  CHECK_THAT(87.96, WithinAbs(L.gamma(), eps));
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_serialization",
                 "[PeriodicStructureAdapter]") {
  StructureAdapterPtr kbise1 = dumpandload(m_kbise);
  CHECK(m_kbise.get() != kbise1.get());
  const double eps = 1.0e-7;
  CHECK(23 == kbise1->countSites());
  CHECK(23.0 == kbise1->totalOccupancy());
  CHECK_THAT(0.0335565, WithinAbs(kbise1->numberDensity(), eps));
  CHECK(string("K1+") == kbise1->siteAtomType(0));
  CHECK(string("Bi3+") == kbise1->siteAtomType(2));
  CHECK(string("Se") == kbise1->siteAtomType(10));
  CHECK(string("Se") == kbise1->siteAtomType(22));
  PeriodicStructureAdapterPtr pkbise =
    std::dynamic_pointer_cast<PeriodicStructureAdapter>(m_kbise);
  PeriodicStructureAdapterPtr pkbise1 =
    std::dynamic_pointer_cast<PeriodicStructureAdapter>(kbise1);
  const Lattice& L = pkbise->getLattice();
  const Lattice& L1 = pkbise1->getLattice();
  CHECK(L.a() == L1.a());
  CHECK(L.b() == L1.b());
  CHECK(L.c() == L1.c());
  CHECK(L.alpha() == L1.alpha());
  CHECK(L.beta() == L1.beta());
  CHECK(L.gamma() == L1.gamma());
}

TEST_CASE_METHOD(TestPeriodicStructureAdapterFixture,
                 "TestPeriodicStructureAdapter::test_comparison",
                 "[PeriodicStructureAdapter]") {
  const PeriodicStructureAdapter& kbise0 =
    static_cast<const PeriodicStructureAdapter&>(*m_kbise);
  PeriodicStructureAdapter kbise1(kbise0);
  PeriodicStructureAdapter kbise2(kbise0);
  CHECK(kbise0 == kbise1);
  CHECK(!(kbise0 != kbise1));
  CHECK(kbise1.countSites());
  kbise1.erase(0);
  CHECK(kbise0 != kbise1);
  CHECK(!(kbise0 == kbise1));
  kbise2.setLatPar(3, 4, 5, 91, 92, 93);
  CHECK(kbise0 != kbise2);
}

//////////////////////////////////////////////////////////////////////////////
// TestPeriodicStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

class TestPeriodicStructureBondGeneratorFixture {
 protected:
  StructureAdapterPtr m_ni;
  BaseBondGeneratorPtr m_nibnds;

 public:
  TestPeriodicStructureBondGeneratorFixture() {
    m_ni = loadTestPeriodicStructure("Ni.stru");
    m_nibnds = m_ni->createBondGenerator();
  }
};  // class TestPeriodicStructureBondGeneratorFixture

TEST_CASE_METHOD(TestPeriodicStructureBondGeneratorFixture,
                 "TestPeriodicStructureBondGenerator::test_typeid",
                 "[PeriodicStructureBondGenerator]") {
  BaseBondGenerator& r_nibnds = *m_nibnds;
  CHECK(typeid(PeriodicStructureBondGenerator) == typeid(r_nibnds));
}

TEST_CASE_METHOD(TestPeriodicStructureBondGeneratorFixture,
                 "TestPeriodicStructureBondGenerator::test_bondCountNickel",
                 "[PeriodicStructureBondGenerator]") {
  m_nibnds->selectAnchorSite(0);
  m_nibnds->setRmin(0);
  m_nibnds->setRmax(1.0);
  CHECK(0 == countBonds(*m_nibnds));
  m_nibnds->selectAnchorSite(3);
  CHECK(0 == countBonds(*m_nibnds));
  m_nibnds->setRmin(-10);
  CHECK(0 == countBonds(*m_nibnds));
  // there are 12 nearest neighbors at 2.49
  m_nibnds->setRmax(3);
  CHECK(12 == countBonds(*m_nibnds));
  m_nibnds->selectAnchorSite(0);
  CHECK(12 == countBonds(*m_nibnds));
  // there are no self neighbors below the cell length of 3.52
  m_nibnds->selectAnchorSite(0);
  m_nibnds->selectSiteRange(0, 1);
  CHECK(0 == countBonds(*m_nibnds));
  // and any other unit cell atom would give 4 neighbors
  m_nibnds->selectAnchorSite(0);
  m_nibnds->selectSiteRange(3, 4);
  CHECK(4 == countBonds(*m_nibnds));
  // there are no bonds between 2.6 and 3.4
  m_nibnds->setRmin(2.6);
  m_nibnds->setRmax(3.4);
  m_nibnds->selectSiteRange(0, 4);
  CHECK(0 == countBonds(*m_nibnds));
  // there are 6 second nearest neighbors at 3.52
  m_nibnds->setRmax(3.6);
  CHECK(6 == countBonds(*m_nibnds));
  // which sums to 18 neighbors within radius 3.6
  m_nibnds->setRmin(0);
  CHECK(18 == countBonds(*m_nibnds));
}

TEST_CASE_METHOD(TestPeriodicStructureBondGeneratorFixture,
                 "TestPeriodicStructureBondGenerator::test_bondCountWurtzite",
                 "[PeriodicStructureBondGenerator]") {
  StructureAdapterPtr stru = loadTestPeriodicStructure("ZnS_wurtzite.stru");
  BaseBondGeneratorPtr bnds = stru->createBondGenerator();
  CHECK(4 == stru->countSites());
  bnds->selectAnchorSite(0);
  bnds->selectSiteRange(0, 4);
  bnds->setRmin(0);
  // there should be no bond below the ZnS distance 2.31
  bnds->setRmax(2.2);
  CHECK(0 == countBonds(*bnds));
  // z-neighbor is slightly more distant than 3 in the lower plane
  bnds->setRmax(2.35);
  CHECK(3 == countBonds(*bnds));
  bnds->setRmax(2.5);
  CHECK(4 == countBonds(*bnds));
  // there are 12 second nearest neighbors at 3.81
  bnds->setRmin(3.7);
  bnds->setRmax(3.82);
  CHECK(12 == countBonds(*bnds));
  // and one more at 3.83
  bnds->setRmax(3.85);
  CHECK(13 == countBonds(*bnds));
  // making the total 17
  bnds->setRmin(0);
  CHECK(17 == countBonds(*bnds));
  // and the same happens for all other sites
  bnds->selectAnchorSite(1);
  CHECK(17 == countBonds(*bnds));
  bnds->selectAnchorSite(2);
  CHECK(17 == countBonds(*bnds));
  bnds->selectAnchorSite(3);
  CHECK(17 == countBonds(*bnds));
}

TEST_CASE_METHOD(TestPeriodicStructureBondGeneratorFixture,
                 "TestPeriodicStructureBondGenerator::test_LiTaO3",
                 "[PeriodicStructureBondGenerator]") {
  const string lithium = "Li1+";
  const string tantalum = "Ta5+";
  const string oxygen = "O2-";
  const double epsu = 1e-5;
  StructureAdapterPtr stru = loadTestPeriodicStructure("LiTaO3.stru");
  CHECK(30 == stru->countSites());
  BaseBondGeneratorPtr bnds = stru->createBondGenerator();
  bnds->selectAnchorSite(0);
  bnds->selectSiteRange(0, 30);
  // there are 3 oxygen neighbors at 2.065
  bnds->setRmax(2.1);
  CHECK(3 == countBonds(*bnds));
  // Li at site 0 is isotropic, oxygens have equal msd-s towards Li
  for (bnds->rewind(); !bnds->finished(); bnds->next()) {
    CHECK(lithium == stru->siteAtomType(bnds->site0()));
    CHECK(oxygen == stru->siteAtomType(bnds->site1()));
    CHECK_THAT(0.00265968, WithinAbs(testmsd0(stru, bnds), epsu));
    CHECK_THAT(0.00710945, WithinAbs(testmsd1(stru, bnds), epsu));
  }
  // there are 3 oxygen neighbors at 2.26
  bnds->setRmin(2.2);
  bnds->setRmax(2.3);
  CHECK(3 == countBonds(*bnds));
  for (bnds->rewind(); !bnds->finished(); bnds->next()) {
    CHECK(oxygen == stru->siteAtomType(bnds->site1()));
    CHECK_THAT(0.00265968, WithinAbs(testmsd0(stru, bnds), epsu));
    CHECK_THAT(0.00824319, WithinAbs(testmsd1(stru, bnds), epsu));
  }
  // finally there are 4 Ta neighbors between 2.8 and 3.1
  bnds->setRmin(2.8);
  bnds->setRmax(3.1);
  CHECK(4 == countBonds(*bnds));
  for (bnds->rewind(); !bnds->finished(); bnds->next()) {
    CHECK_THAT(0.00265968, WithinAbs(testmsd0(stru, bnds), epsu));
    CHECK(tantalum == stru->siteAtomType(bnds->site1()));
    R3::Vector r01xy = bnds->r01();
    r01xy[2] = 0.0;
    // for the tantalum above Li the msd equals U33
    if (R3::norm(r01xy) < 0.1) {
      CHECK_THAT(0.00356, WithinAbs(testmsd1(stru, bnds), epsu));
    }
    // other 3 tantalums are related by triple axis and
    // have the same msd towards the central Li
    else {
      CHECK_THAT(0.00486942, WithinAbs(testmsd1(stru, bnds), epsu));
    }
  }
}

}  // namespace srreal
}  // namespace diffpy

// End of file
