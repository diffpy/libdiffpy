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
 * TestObjCrystStructureAdapter -- unit tests for an adapter
 *     to the ObjCryst::Crystal class
 *
 * TestObjCrystStructureBondGenerator -- unit tests for bond generator
 *
 *****************************************************************************/

#include <typeinfo>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include <diffpy/srreal/StructureDifference.hpp>
#include "serialization_helpers.hpp"
#include "objcryst_helpers.hpp"
#include <ObjCryst/ObjCryst/Crystal.h>
#include <ObjCryst/ObjCryst/Molecule.h>
#include <ObjCryst/ObjCryst/ScatteringPower.h>

using Catch::Matchers::WithinAbs;

using namespace std;
using namespace diffpy::srreal;
using ObjCryst::Crystal;
using ObjCryst::Molecule;
using ObjCryst::ScatteringPowerAtom;

// Local Helpers -------------------------------------------------------------

namespace {

int countBonds(diffpy::srreal::BaseBondGenerator& bnds) {
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

double diffdegree(StructureAdapterPtr stru0, StructureAdapterPtr stru1) {
  StructureDifference sd = stru0->diff(stru1);
  double tot = 0.0;
  if (stru0) tot += stru0->countSites();
  if (stru1) tot += stru1->countSites();
  double rv = (sd.pop0.size() + sd.add1.size()) / (tot ? tot : 1.0);
  return rv;
}

// Defined below
Molecule* makeC60Molecule();

}  // namespace

//////////////////////////////////////////////////////////////////////////////
// TestObjCrystStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestObjCrystStructureAdapterFixture {
 protected:
  StructureAdapterPtr m_ni =
    createStructureAdapter(cachedTestCrystal("Ni.cif"));
  StructureAdapterPtr m_kbise =
    createStructureAdapter(cachedTestCrystal("alpha_K2Bi8Se13.cif"));
  StructureAdapterPtr m_catio3 =
    createStructureAdapter(cachedTestCrystal("CaTiO3.cif"));
};  // class TestObjCrystStructureAdapterFixture

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_typeid",
                 "[ObjCrystStructureAdapter]") {
  StructureAdapter& r_ni = *m_ni;
  CHECK(typeid(CrystalStructureAdapter) == typeid(r_ni));
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_countSites",
                 "[ObjCrystStructureAdapter]") {
  CHECK(1 == m_ni->countSites());
  CHECK(12 == m_kbise->countSites());
  CHECK(4 == m_catio3->countSites());
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_totalOccupancy",
                 "[ObjCrystStructureAdapter]") {
  CHECK(4.0 == m_ni->totalOccupancy());
  CHECK(23.0 == m_kbise->totalOccupancy());
  CHECK(20.0 == m_catio3->totalOccupancy());
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_numberDensity",
                 "[ObjCrystStructureAdapter]") {
  const double eps = 1.0e-7;
  CHECK_THAT(0.0914114, WithinAbs(m_ni->numberDensity(), eps));
  CHECK_THAT(0.0335565, WithinAbs(m_kbise->numberDensity(), eps));
  CHECK_THAT(0.0894566, WithinAbs(m_catio3->numberDensity(), eps));
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_siteCartesianPosition",
                 "[ObjCrystStructureAdapter]") {
  const double eps = 1.0e-5;
  R3::Vector rCa = m_catio3->siteCartesianPosition(1);
  CHECK_THAT(0.03486, WithinAbs(rCa[0], eps));
  CHECK_THAT(0.19366, WithinAbs(rCa[1], eps));
  CHECK_THAT(1.90975, WithinAbs(rCa[2], eps));
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_siteAnisotropy",
                 "[ObjCrystStructureAdapter]") {
  for (int i = 0; i < m_ni->countSites(); ++i) {
    CHECK(false == m_ni->siteAnisotropy(i));
  }
  for (int i = 0; i < m_catio3->countSites(); ++i) {
    CHECK(true == m_catio3->siteAnisotropy(i));
  }
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_siteCartesianUij",
                 "[ObjCrystStructureAdapter]") {
  // nickel should have all Uij equal zero.
  const double* puij = &(m_ni->siteCartesianUij(0).data()[0]);
  for (int i = 0; i < 9; ++i) {
    CHECK(0.0 == puij[i]);
  }
  // check CaTiO3 values
  const R3::Matrix& UO2 = m_catio3->siteCartesianUij(3);
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

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_siteAtomType",
                 "[ObjCrystStructureAdapter]") {
  CHECK(string("Ni") == m_ni->siteAtomType(0));
  CHECK(string("K1+") == m_kbise->siteAtomType(0));
  CHECK(string("Bi3+") == m_kbise->siteAtomType(1));
  CHECK(string("Se") == m_kbise->siteAtomType(5));
  CHECK(string("Se") == m_kbise->siteAtomType(11));
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_diff",
                 "[ObjCrystStructureAdapter]") {
  StructureAdapterPtr acto = m_catio3;
  CrystalStructureAdapterPtr ccto(
    std::dynamic_pointer_cast<CrystalStructureAdapter>(acto));
  StructureAdapterPtr acto2 = acto->clone();
  CrystalStructureAdapterPtr ccto2(
    std::dynamic_pointer_cast<CrystalStructureAdapter>(acto2));
  CHECK(0.0 == diffdegree(acto, acto));
  CHECK(1.0 == diffdegree(acto, StructureAdapterPtr()));
  CHECK(0.0 == diffdegree(acto, acto2));
  CHECK(0.0 == diffdegree(acto, ccto2));
  const double symprec = ccto2->getSymmetryPrecision();
  ccto2->setSymmetryPrecision(0.010723);
  CHECK(1.0 == diffdegree(acto, acto2));
  ccto2->setSymmetryPrecision(symprec);
  CHECK(0.0 == diffdegree(acto, acto2));
  ccto2->clearSymOps();
  CHECK(1.0 == diffdegree(acto, acto2));
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_getLattice",
                 "[ObjCrystStructureAdapter]") {
  CrystalStructureAdapter* pkbise =
    dynamic_cast<CrystalStructureAdapter*>(m_kbise.get());
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

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_getEquivalentAtoms",
                 "[ObjCrystStructureAdapter]") {
  diffpy::mathutils::EpsilonEqual allclose;
  unique_ptr<Crystal> crst(loadTestCrystal("NH4Br.cif"));
  StructureAdapterPtr adpt = createStructureAdapter(*crst);
  CrystalStructureAdapterPtr cadpt =
    std::dynamic_pointer_cast<CrystalStructureAdapter>(adpt);
  CHECK(2 == cadpt->countSites());
  CHECK(4 == cadpt->totalOccupancy());
  CrystalStructureAdapter::AtomVector eq0, eq1;
  eq0 = cadpt->getEquivalentAtoms(0);
  CHECK(2 == eq0.size());
  cadpt->toFractional(eq0[0]);
  cadpt->toFractional(eq0[1]);
  CHECK(allclose(R3::Vector(0, 0.5, 0.03), eq0[0].xyz_cartn));
  CHECK(allclose(R3::Vector(0.5, 0, 0.97), eq0[1].xyz_cartn));
  eq1 = cadpt->getEquivalentAtoms(1);
  CHECK(2 == eq1.size());
  cadpt->toFractional(eq1[0]);
  cadpt->toFractional(eq1[1]);
  CHECK(allclose(R3::Vector(0, 0, 0.5), eq1[0].xyz_cartn));
  CHECK(allclose(R3::Vector(0.5, 0.5, 0.5), eq1[1].xyz_cartn));
}

TEST_CASE_METHOD(TestObjCrystStructureAdapterFixture,
                 "TestObjCrystStructureAdapter::test_serialization",
                 "[ObjCrystStructureAdapter]") {
  stringstream storage(ios::in | ios::out | ios::binary);
  diffpy::serialization::oarchive oa(storage, ios::binary);
  oa << m_kbise;
  diffpy::serialization::iarchive ia(storage, ios::binary);
  StructureAdapterPtr kbise1;
  CHECK(!kbise1.get());
  ia >> kbise1;
  CHECK(m_kbise.get() != kbise1.get());
  const double eps = 1.0e-7;
  CHECK(12 == kbise1->countSites());
  CHECK(23.0 == kbise1->totalOccupancy());
  CHECK_THAT(0.0335565, WithinAbs(kbise1->numberDensity(), eps));
  CHECK(string("K1+") == kbise1->siteAtomType(0));
  CHECK(string("Bi3+") == kbise1->siteAtomType(1));
  CHECK(string("Se") == kbise1->siteAtomType(5));
  CHECK(string("Se") == kbise1->siteAtomType(11));
  CrystalStructureAdapter* pkbise =
    dynamic_cast<CrystalStructureAdapter*>(m_kbise.get());
  CrystalStructureAdapter* pkbise1 =
    dynamic_cast<CrystalStructureAdapter*>(kbise1.get());
  const Lattice& L = pkbise->getLattice();
  const Lattice& L1 = pkbise1->getLattice();
  CHECK(L.a() == L1.a());
  CHECK(L.b() == L1.b());
  CHECK(L.c() == L1.c());
  CHECK(L.alpha() == L1.alpha());
  CHECK(L.beta() == L1.beta());
  CHECK(L.gamma() == L1.gamma());
}

//////////////////////////////////////////////////////////////////////////////
// TestObjCrystStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

class TestObjCrystStructureBondGeneratorFixture {
 protected:
  StructureAdapterPtr m_ni =
    createStructureAdapter(cachedTestCrystal("Ni.cif"));
  BaseBondGeneratorPtr m_nibnds = m_ni->createBondGenerator();
};  // class TestObjCrystStructureBondGeneratorFixture

TEST_CASE_METHOD(TestObjCrystStructureBondGeneratorFixture,
                 "TestObjCrystStructureBondGenerator::test_typeid",
                 "[ObjCrystStructureBondGenerator]") {
  BaseBondGenerator& r_nibnds = *m_nibnds;
  CHECK(typeid(CrystalStructureBondGenerator) == typeid(r_nibnds));
}

TEST_CASE_METHOD(TestObjCrystStructureBondGeneratorFixture,
                 "TestObjCrystStructureBondGenerator::test_bondCountNickel",
                 "[ObjCrystStructureBondGenerator]") {
  m_nibnds->selectAnchorSite(0);
  m_nibnds->setRmin(0);
  m_nibnds->setRmax(1.0);
  CHECK(0 == countBonds(*m_nibnds));
  m_nibnds->setRmin(-10);
  CHECK(0 == countBonds(*m_nibnds));
  // there are 12 nearest neighbors at 2.49
  m_nibnds->setRmax(3);
  CHECK(12 == countBonds(*m_nibnds));
  // there are no bonds between 2.6 and 3.4
  m_nibnds->setRmin(2.6);
  m_nibnds->setRmax(3.4);
  CHECK(0 == countBonds(*m_nibnds));
  // there are 6 second nearest neighbors at 3.52
  m_nibnds->setRmax(3.6);
  CHECK(6 == countBonds(*m_nibnds));
  // which sums to 18 neighbors within radius 3.6
  m_nibnds->setRmin(0);
  CHECK(18 == countBonds(*m_nibnds));
}

TEST_CASE_METHOD(TestObjCrystStructureBondGeneratorFixture,
                 "TestObjCrystStructureBondGenerator::test_bondCountWurtzite",
                 "[ObjCrystStructureBondGenerator]") {
  unique_ptr<Crystal> wurtzite(loadTestCrystal("ZnS_wurtzite.cif"));
  StructureAdapterPtr stru = createStructureAdapter(*wurtzite);
  BaseBondGeneratorPtr bnds = stru->createBondGenerator();
  CHECK(2 == stru->countSites());
  bnds->selectAnchorSite(0);
  bnds->selectSiteRange(0, 2);
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
}

TEST_CASE_METHOD(TestObjCrystStructureBondGeneratorFixture,
                 "TestObjCrystStructureBondGenerator::test_LiTaO3",
                 "[ObjCrystStructureBondGenerator]") {
  const string lithium = "Li1+";
  const string tantalum = "Ta5+";
  const string oxygen = "O2-";
  const double epsu = 1e-5;
  unique_ptr<Crystal> litao3(loadTestCrystal("LiTaO3.cif"));
  StructureAdapterPtr stru = createStructureAdapter(*litao3);
  CHECK(3 == stru->countSites());
  BaseBondGeneratorPtr bnds = stru->createBondGenerator();
  bnds->selectAnchorSite(0);
  bnds->selectSiteRange(0, 3);
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

//////////////////////////////////////////////////////////////////////////////
// TestObjCrystMoleculeAdapter
//////////////////////////////////////////////////////////////////////////////

class TestObjCrystMoleculeAdapterFixture {
 protected:
  Molecule* mmol_c60;
  StructureAdapterPtr m_c60;

 public:
  TestObjCrystMoleculeAdapterFixture() {
    mmol_c60 = makeC60Molecule();
    m_c60 = createStructureAdapter(*mmol_c60);
  }

  ~TestObjCrystMoleculeAdapterFixture() {
    Crystal* crst = &(mmol_c60->GetCrystal());
    delete crst;
  }
};  // class TestObjCrystMoleculeAdapterFixture

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_typeid",
                 "[ObjCrystMoleculeAdapter]") {
  StructureAdapter& r_c60 = *m_c60;
  CHECK(typeid(AtomicStructureAdapter) == typeid(r_c60));
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_countSites",
                 "[ObjCrystMoleculeAdapter]") {
  CHECK(60 == m_c60->countSites());
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_totalOccupancy",
                 "[ObjCrystMoleculeAdapter]") {
  CHECK(60 == m_c60->totalOccupancy());
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_numberDensity",
                 "[ObjCrystMoleculeAdapter]") {
  const double eps = 1.0e-7;
  CHECK_THAT(0.0, WithinAbs(m_c60->numberDensity(), eps));
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_siteCartesianPosition",
                 "[ObjCrystMoleculeAdapter]") {
  const double eps = 1.0e-5;
  R3::Vector rC0 = m_c60->siteCartesianPosition(0);
  CHECK_THAT(3.45127, WithinAbs(rC0[0], eps));
  CHECK_THAT(0.68500, WithinAbs(rC0[1], eps));
  CHECK_THAT(0.00000, WithinAbs(rC0[2], eps));
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_siteAnisotropy",
                 "[ObjCrystMoleculeAdapter]") {
  for (int i = 0; i < m_c60->countSites(); ++i) {
    CHECK(false == m_c60->siteAnisotropy(i));
  }
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_siteCartesianUij",
                 "[ObjCrystMoleculeAdapter]") {
  double BtoU = 1.0 / (8 * M_PI * M_PI);
  const R3::Matrix Uij = m_c60->siteCartesianUij(0);
  CHECK(0.0 == Uij(0, 1));
  CHECK(0.0 == Uij(0, 2));
  CHECK(0.0 == Uij(1, 0));
  CHECK(0.0 == Uij(1, 2));
  CHECK(0.0 == Uij(2, 0));
  CHECK(0.0 == Uij(2, 1));
  CHECK(BtoU == Uij(0, 0));
  CHECK(BtoU == Uij(1, 1));
  CHECK(BtoU == Uij(2, 2));
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_siteAtomType",
                 "[ObjCrystMoleculeAdapter]") {
  for (int i = 0; i < m_c60->countSites(); ++i) {
    CHECK(string("C") == m_c60->siteAtomType(i));
  }
}

TEST_CASE_METHOD(TestObjCrystMoleculeAdapterFixture,
                 "TestObjCrystMoleculeAdapter::test_serialization",
                 "[ObjCrystMoleculeAdapter]") {
  StructureAdapterPtr c60a;
  c60a = dumpandload(m_c60);
  CHECK(m_c60.get() != c60a.get());
  CHECK(60 == c60a->countSites());
  CHECK(60.0 == c60a->totalOccupancy());
  CHECK(0.0 == c60a->numberDensity());
  R3::Vector dxyz(1, 2, 3);
  dxyz = m_c60->siteCartesianPosition(0) - c60a->siteCartesianPosition(0);
  CHECK(0.0 == R3::norm(dxyz));
  dxyz = m_c60->siteCartesianPosition(37) - c60a->siteCartesianPosition(37);
  CHECK(0.0 == R3::norm(dxyz));
  for (int i = 0; i < c60a->countSites(); ++i) {
    CHECK(string("C") == c60a->siteAtomType(i));
  }
}

//////////////////////////////////////////////////////////////////////////////
// TestObjCrystMoleculeBondGenerator
//////////////////////////////////////////////////////////////////////////////

class TestObjCrystMoleculeBondGeneratorFixture {
 protected:
  Molecule* mmol_c60;
  StructureAdapterPtr m_c60;
  BaseBondGeneratorPtr m_c60bnds;

 public:
  TestObjCrystMoleculeBondGeneratorFixture() {
    mmol_c60 = makeC60Molecule();
    m_c60 = createStructureAdapter(*mmol_c60);
    m_c60bnds = m_c60->createBondGenerator();
  }

  ~TestObjCrystMoleculeBondGeneratorFixture() {
    Crystal* crst = &(mmol_c60->GetCrystal());
    delete crst;
  }
};  // class TestObjCrystMoleculeBondGeneratorFixture

TEST_CASE_METHOD(TestObjCrystMoleculeBondGeneratorFixture,
                 "TestObjCrystMoleculeBondGenerator::test_typeid",
                 "[ObjCrystMoleculeBondGenerator]") {
  BaseBondGenerator& r_c60bnds = *m_c60bnds;
  CHECK(typeid(BaseBondGenerator) == typeid(r_c60bnds));
}

TEST_CASE_METHOD(TestObjCrystMoleculeBondGeneratorFixture,
                 "TestObjCrystMoleculeBondGenerator::test_bondCountC60",
                 "[ObjCrystMoleculeBondGenerator]") {
  m_c60bnds->selectAnchorSite(0);
  m_c60bnds->setRmin(0);
  m_c60bnds->setRmax(1.0);
  CHECK(0 == countBonds(*m_c60bnds));
  m_c60bnds->setRmin(-10);
  CHECK(0 == countBonds(*m_c60bnds));
  // there are 3 nearest neighbors at 1.44
  m_c60bnds->setRmax(1.5);
  CHECK(3 == countBonds(*m_c60bnds));
  // There are 59 neighbors total
  m_c60bnds->setRmax(8);
  CHECK(59 == countBonds(*m_c60bnds));
}

namespace {

Molecule* makeC60Molecule() {
  Crystal* cryst = new Crystal(1, 1, 1, "P1");
  Molecule* mol = new Molecule(*cryst, "c60");
  cryst->AddScatterer(mol);

  // Populate the molecule
  ScatteringPowerAtom* sp = new ScatteringPowerAtom("C", "C");
  cryst->AddScatteringPower(sp);
  mol->AddAtom(3.451266498, 0.685000000, 0.000000000, sp, "C0");
  mol->AddAtom(3.451266498, -0.685000000, 0.000000000, sp, "C1");
  mol->AddAtom(-3.451266498, 0.685000000, 0.000000000, sp, "C2");
  mol->AddAtom(-3.451266498, -0.685000000, 0.000000000, sp, "C3");
  mol->AddAtom(0.685000000, 0.000000000, 3.451266498, sp, "C4");
  mol->AddAtom(-0.685000000, 0.000000000, 3.451266498, sp, "C5");
  mol->AddAtom(0.685000000, 0.000000000, -3.451266498, sp, "C6");
  mol->AddAtom(-0.685000000, 0.000000000, -3.451266498, sp, "C7");
  mol->AddAtom(0.000000000, 3.451266498, 0.685000000, sp, "C8");
  mol->AddAtom(0.000000000, 3.451266498, -0.685000000, sp, "C9");
  mol->AddAtom(0.000000000, -3.451266498, 0.685000000, sp, "C10");
  mol->AddAtom(0.000000000, -3.451266498, -0.685000000, sp, "C11");
  mol->AddAtom(3.003809890, 1.409000000, 1.171456608, sp, "C12");
  mol->AddAtom(3.003809890, 1.409000000, -1.171456608, sp, "C13");
  mol->AddAtom(3.003809890, -1.409000000, 1.171456608, sp, "C14");
  mol->AddAtom(3.003809890, -1.409000000, -1.171456608, sp, "C15");
  mol->AddAtom(-3.003809890, 1.409000000, 1.171456608, sp, "C16");
  mol->AddAtom(-3.003809890, 1.409000000, -1.171456608, sp, "C17");
  mol->AddAtom(-3.003809890, -1.409000000, 1.171456608, sp, "C18");
  mol->AddAtom(-3.003809890, -1.409000000, -1.171456608, sp, "C19");
  mol->AddAtom(1.409000000, 1.171456608, 3.003809890, sp, "C20");
  mol->AddAtom(1.409000000, -1.171456608, 3.003809890, sp, "C21");
  mol->AddAtom(-1.409000000, 1.171456608, 3.003809890, sp, "C22");
  mol->AddAtom(-1.409000000, -1.171456608, 3.003809890, sp, "C23");
  mol->AddAtom(1.409000000, 1.171456608, -3.003809890, sp, "C24");
  mol->AddAtom(1.409000000, -1.171456608, -3.003809890, sp, "C25");
  mol->AddAtom(-1.409000000, 1.171456608, -3.003809890, sp, "C26");
  mol->AddAtom(-1.409000000, -1.171456608, -3.003809890, sp, "C27");
  mol->AddAtom(1.171456608, 3.003809890, 1.409000000, sp, "C28");
  mol->AddAtom(-1.171456608, 3.003809890, 1.409000000, sp, "C29");
  mol->AddAtom(1.171456608, 3.003809890, -1.409000000, sp, "C30");
  mol->AddAtom(-1.171456608, 3.003809890, -1.409000000, sp, "C31");
  mol->AddAtom(1.171456608, -3.003809890, 1.409000000, sp, "C32");
  mol->AddAtom(-1.171456608, -3.003809890, 1.409000000, sp, "C33");
  mol->AddAtom(1.171456608, -3.003809890, -1.409000000, sp, "C34");
  mol->AddAtom(-1.171456608, -3.003809890, -1.409000000, sp, "C35");
  mol->AddAtom(2.580456608, 0.724000000, 2.279809890, sp, "C36");
  mol->AddAtom(2.580456608, 0.724000000, -2.279809890, sp, "C37");
  mol->AddAtom(2.580456608, -0.724000000, 2.279809890, sp, "C38");
  mol->AddAtom(2.580456608, -0.724000000, -2.279809890, sp, "C39");
  mol->AddAtom(-2.580456608, 0.724000000, 2.279809890, sp, "C40");
  mol->AddAtom(-2.580456608, 0.724000000, -2.279809890, sp, "C41");
  mol->AddAtom(-2.580456608, -0.724000000, 2.279809890, sp, "C42");
  mol->AddAtom(-2.580456608, -0.724000000, -2.279809890, sp, "C43");
  mol->AddAtom(0.724000000, 2.279809890, 2.580456608, sp, "C44");
  mol->AddAtom(0.724000000, -2.279809890, 2.580456608, sp, "C45");
  mol->AddAtom(-0.724000000, 2.279809890, 2.580456608, sp, "C46");
  mol->AddAtom(-0.724000000, -2.279809890, 2.580456608, sp, "C47");
  mol->AddAtom(0.724000000, 2.279809890, -2.580456608, sp, "C48");
  mol->AddAtom(0.724000000, -2.279809890, -2.580456608, sp, "C49");
  mol->AddAtom(-0.724000000, 2.279809890, -2.580456608, sp, "C50");
  mol->AddAtom(-0.724000000, -2.279809890, -2.580456608, sp, "C51");
  mol->AddAtom(2.279809890, 2.580456608, 0.724000000, sp, "C52");
  mol->AddAtom(-2.279809890, 2.580456608, 0.724000000, sp, "C53");
  mol->AddAtom(2.279809890, 2.580456608, -0.724000000, sp, "C54");
  mol->AddAtom(-2.279809890, 2.580456608, -0.724000000, sp, "C55");
  mol->AddAtom(2.279809890, -2.580456608, 0.724000000, sp, "C56");
  mol->AddAtom(-2.279809890, -2.580456608, 0.724000000, sp, "C57");
  mol->AddAtom(2.279809890, -2.580456608, -0.724000000, sp, "C58");
  mol->AddAtom(-2.279809890, -2.580456608, -0.724000000, sp, "C59");

  return mol;
}

}  // namespace

// End of file
