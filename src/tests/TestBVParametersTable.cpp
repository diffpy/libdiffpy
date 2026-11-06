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
 * TestBVParametersTable -- unit tests for BVParametersTable class
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/BVParametersTable.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

class TestBVParametersTableFixture {
 protected:
  BVParametersTablePtr mbvtb;

 public:
  TestBVParametersTableFixture() { mbvtb.reset(new BVParametersTable); }
};  // class TestBVParametersTableFixture

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_copy", "[BVParametersTable]") {
  BVParam mynacl("Na", 1, "Cl", -1);
  mynacl.mRo = 1.23;
  mynacl.mB = 0.377;
  mbvtb->setCustom(mynacl);
  CHECK(1.23 == mbvtb->lookup("Na", 1, "Cl", -1).mRo);
  CHECK(0.377 == mbvtb->lookup("Na", 1, "Cl", -1).mB);
  BVParametersTable bvcp(*mbvtb);
  CHECK(1.23 == bvcp.lookup("Na", 1, "Cl", -1).mRo);
  CHECK(0.377 == bvcp.lookup("Na", 1, "Cl", -1).mB);
  bvcp.resetAll();
  CHECK(2.15 == bvcp.lookup("Na", 1, "Cl", -1).mRo);
  CHECK(0.37 == bvcp.lookup("Na", 1, "Cl", -1).mB);
  CHECK(1.23 == mbvtb->lookup("Na", 1, "Cl", -1).mRo);
  CHECK(0.377 == mbvtb->lookup("Na", 1, "Cl", -1).mB);
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_none", "[BVParametersTable]") {
  const BVParam& bpnone = BVParametersTable::none();
  CHECK(&bpnone == &(mbvtb->none()));
  CHECK(0.0 == bpnone.bondvalence(0));
  CHECK(0.0 == bpnone.bondvalence(3));
  CHECK(0.0 == bpnone.bondvalence(10));
  CHECK(0.0 == bpnone.bondvalenceToDistance(0.5));
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_atomvalence",
                 "[BVParametersTable]") {
  CHECK(0 == mbvtb->getAtomValence("O"));
  mbvtb->setAtomValence("O", -2);
  CHECK(-2 == mbvtb->getAtomValence("O"));
  mbvtb->resetAtomValences();
  CHECK(0 == mbvtb->getAtomValence("O"));
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_lookup", "[BVParametersTable]") {
  BVParam bp = mbvtb->lookup("Xx", 0, "Yy", 3);
  CHECK(bp.matom0.empty());
  CHECK(bp.matom1.empty());
  CHECK(0 == bp.mvalence0);
  CHECK(0 == bp.mvalence1);
  CHECK(0.0 == bp.mRo);
  CHECK(0.0 == bp.mB);
  CHECK(&mbvtb->lookup("Xx", 0, "Yy", 3) ==
        &mbvtb->lookup("Na", 33, "Cl", -11));
  BVParam bk("Na", 1, "Cl", -1);
  CHECK(&mbvtb->lookup(bk) == &mbvtb->lookup("Cl", -1, "Na", 1));
  CHECK(2.15 == mbvtb->lookup(bk).mRo);
  CHECK(0.37 == mbvtb->lookup(bk).mB);
  CHECK("b" == mbvtb->lookup(bk).mref_id);
  const BVParam& bnacl0 = mbvtb->lookup("Na", 1, "Cl", -1);
  CHECK(&bnacl0 == &mbvtb->lookup("Na+", "Cl-"));
  CHECK(&bnacl0 == &mbvtb->lookup("Cl1-", "Na1+"));
  const BVParam& bpnone = BVParametersTable::none();
  CHECK(&bpnone == &mbvtb->lookup("Na", "Cl"));
  mbvtb->setAtomValence("Na", +1);
  mbvtb->setAtomValence("Cl", -1);
  CHECK(&bnacl0 == &mbvtb->lookup("Na", "Cl"));
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_setCustom",
                 "[BVParametersTable]") {
  BVParam mymgo("Mg", 2, "O", -2);
  mymgo.mRo = 2.34;
  mymgo.mB = 0.345;
  mbvtb->setCustom(mymgo);
  CHECK(2.34 == mbvtb->lookup("Mg", 2, "O", -2).mRo);
  CHECK(0.345 == mbvtb->lookup(mymgo).mB);
  mbvtb->resetCustom(mymgo);
  CHECK(2.34 != mbvtb->lookup("Mg", 2, "O", -2).mRo);
  CHECK(0.345 != mbvtb->lookup(mymgo).mB);
  BVParam zrh = mbvtb->lookup("Zr", 4, "H", -1);
  CHECK(1.79 == zrh.mRo);
  CHECK(0.37 == zrh.mB);
  CHECK(string("b") == zrh.mref_id);
  mbvtb->setCustom("H", -1, "Zr", 4, 1.791, 0.371, "check");
  zrh = mbvtb->lookup("Zr", 4, "H", -1);
  CHECK(1.791 == zrh.mRo);
  CHECK(0.371 == zrh.mB);
  CHECK(string("check") == zrh.mref_id);
  mbvtb->setCustom("H", -1, "Zr", 4, 1.791, 0.444, "mate");
  zrh = mbvtb->lookup("Zr", 4, "H", -1);
  CHECK(0.444 == zrh.mB);
  CHECK(string("mate") == zrh.mref_id);
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_resetAll",
                 "[BVParametersTable]") {
  BVParam mynacl("Cl", -1, "Na", 1, 2.345, 0.44, "pj1");
  BVParam mymgo("O", -2, "Mg", 2, 3.456, 0.55, "pj2");
  mbvtb->setCustom(mynacl);
  mbvtb->setCustom(mymgo);
  CHECK(string("pj1") == mbvtb->lookup(mynacl).mref_id);
  CHECK(string("pj2") == mbvtb->lookup(mymgo).mref_id);
  mbvtb->resetAll();
  CHECK(string("pj1") != mbvtb->lookup(mynacl).mref_id);
  CHECK(string("pj2") != mbvtb->lookup(mymgo).mref_id);
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_getAllCustom",
                 "[BVParametersTable]") {
  CHECK(mbvtb->getAllCustom().empty());
  BVParam mymgo("Mg", 2, "O", -2);
  mbvtb->setCustom(mymgo);
  CHECK(1u == mbvtb->getAllCustom().size());
  const BVParam& bp = mbvtb->lookup("O2-", "Mg2+");
  BVParametersTable::SetOfBVParam::const_iterator ii;
  ii = mbvtb->getAllCustom().find(mymgo);
  CHECK(ii != mbvtb->getAllCustom().end());
  CHECK(&mymgo != &(*ii));
  CHECK(&bp == &(*ii));
  CHECK(mymgo == *ii);
  // erase the only custom parameter
  mbvtb->resetCustom(*ii);
  CHECK(mbvtb->getAllCustom().empty());
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_getAll", "[BVParametersTable]") {
  BVParametersTable::SetOfBVParam allpars0, allpars1;
  allpars0 = mbvtb->getAll();
  CHECK(!allpars0.empty());
  size_t cnt0 = allpars0.size();
  BVParam mymgo("Mg", 2, "O", -2);
  mbvtb->setCustom(mymgo);
  CHECK(cnt0 == mbvtb->getAll().size());
  mbvtb->setCustom(BVParam("Mg", 7, "O", -3));
  allpars1 = mbvtb->getAll();
  size_t cnt1 = allpars1.size();
  CHECK(cnt0 + 1 == cnt1);
  mbvtb->resetAll();
  CHECK(cnt0 == mbvtb->getAll().size());
  mbvtb->setCustom(BVParam("Mg", 2, "O", -2, 77, 88));
  allpars1 = mbvtb->getAll();
  const BVParam& mgo = *allpars1.find(BVParam("Mg", 2, "O", -2));
  CHECK(77.0 == mgo.mRo);
  CHECK(88.0 == mgo.mB);
}

TEST_CASE_METHOD(TestBVParametersTableFixture,
                 "TestBVParametersTable::test_serialization",
                 "[BVParametersTable]") {
  BVParam mynacl("Cl", -1, "Na", 1, 2.345, 0.44, "pj1");
  BVParam mymgo("O", -2, "Mg", 2, 3.456, 0.55, "pj2");
  mbvtb->setCustom(mynacl);
  mbvtb->setCustom(mymgo);
  // check serialization of shared pointers
  BVParametersTablePtr bvtb1 = dumpandload(mbvtb);
  CHECK(mbvtb.get() != bvtb1.get());
  CHECK(2.345 == bvtb1->lookup("Cl", -1, "Na", 1).mRo);
  CHECK(0.44 == bvtb1->lookup("Cl", -1, "Na", 1).mB);
  CHECK(string("pj1") == bvtb1->lookup("Cl", -1, "Na", 1).mref_id);
  CHECK(3.456 == bvtb1->lookup("O", -2, "Mg", 2).mRo);
  CHECK(0.55 == bvtb1->lookup("O", -2, "Mg", 2).mB);
  CHECK(string("pj2") == bvtb1->lookup("O", -2, "Mg", 2).mref_id);
  // check serialization of instances
  BVParametersTable tb2 = dumpandload(*mbvtb);
  CHECK(2u == tb2.getAllCustom().size());
  CHECK(mynacl == tb2.lookup("Cl-", "Na+"));
  // check serialization of customized valences
  tb2.setAtomValence("Na", 1);
  tb2.setAtomValence("Cl", -1);
  BVParametersTable tb3 = dumpandload(tb2);
  CHECK(mynacl == tb3.lookup("Cl", "Na"));
  const BVParam& bpnone = BVParametersTable::none();
  CHECK(bpnone == mbvtb->lookup("Cl", "Na"));
}

// End of file
