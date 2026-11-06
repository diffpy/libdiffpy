/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2010 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * TestNoMetaStructureAdapter -- unit tests for an adapter that
 *     avoids metadata forwarding from the diffpy.Structure classes.
 *
 *****************************************************************************/

#include <typeinfo>
#include <sstream>
#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/NoMetaStructureAdapter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/serialization.hpp>
#include "test_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

//////////////////////////////////////////////////////////////////////////////
// TestNoMetaStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestNoMetaStructureAdapterFixture {
 protected:
  StructureAdapterPtr m_pswt;

 public:
  TestNoMetaStructureAdapterFixture() {
    m_pswt = loadTestPeriodicStructure("PbScW25TiO3.stru");
  }
};  // class TestNoMetaStructureAdapterFixture

TEST_CASE_METHOD(
  TestNoMetaStructureAdapterFixture,
  "TestNoMetaStructureAdapter::test_PDFCalcPreservesConfiguration",
  "[NoMetaStructureAdapter]") {
  PDFCalculator pdfc;
  pdfc.setRstep(0.1);
  pdfc.setRmax(5.0);
  pdfc.setDoubleAttr("scale", 2.5);
  pdfc.setStructure(m_pswt);
  CHECK(2.5 == pdfc.getDoubleAttr("scale"));
  pdfc.eval(m_pswt);
  CHECK(2.5 == pdfc.getDoubleAttr("scale"));
}

TEST_CASE_METHOD(TestNoMetaStructureAdapterFixture,
                 "TestNoMetaStructureAdapter::test_PDFCalcNoMeta",
                 "[NoMetaStructureAdapter]") {
  PDFCalculator pdfc;
  pdfc.setRstep(0.1);
  pdfc.setRmax(5.0);
  CHECK(1.0 == pdfc.getDoubleAttr("scale"));
  pdfc.eval(nometa(m_pswt));
  CHECK(1.0 == pdfc.getDoubleAttr("scale"));
}

TEST_CASE_METHOD(TestNoMetaStructureAdapterFixture,
                 "TestNoMetaStructureAdapter::test_NoMetaTwice",
                 "[NoMetaStructureAdapter]") {
  StructureAdapterPtr adpt0 = nometa(m_pswt);
  StructureAdapterPtr adpt1 = nometa(adpt0);
  CHECK(m_pswt == adpt0);
  CHECK(adpt0 == adpt1);
  CHECK(emptyStructureAdapter() == nometa(StructureAdapterPtr()));
  StructureAdapterPtr proxy(new NoMetaStructureAdapter(m_pswt));
  CHECK(proxy == nometa(proxy));
}

TEST_CASE_METHOD(TestNoMetaStructureAdapterFixture,
                 "TestNoMetaStructureAdapter::test_serialization",
                 "[NoMetaStructureAdapter]") {
  using namespace diffpy::serialization;
  stringstream storage(ios::in | ios::out | ios::binary);
  oarchive oa(storage, ios::binary);
  StructureAdapterPtr pswtbare0(new NoMetaStructureAdapter(m_pswt));
  oa << pswtbare0;
  iarchive ia(storage, ios::binary);
  StructureAdapterPtr pswtbare1;
  CHECK(!pswtbare1.get());
  ia >> pswtbare1;
  CHECK(m_pswt.get() != pswtbare1.get());
  CHECK(56 == pswtbare1->countSites());
  CHECK(string("Pb") == pswtbare1->siteAtomType(0));
  CHECK(string("Ti") == pswtbare1->siteAtomType(55));
  StructureAdapter& r_pswtbare1 = *pswtbare1;
  CHECK(typeid(NoMetaStructureAdapter) == typeid(r_pswtbare1));
}

// End of file
