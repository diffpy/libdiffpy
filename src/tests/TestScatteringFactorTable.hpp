/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestScatteringFactorTable -- unit tests for implementations
*     of the ScatteringFactorTable class
*
* $Id$
*
*****************************************************************************/

#include <typeinfo>
#include <stdexcept>
#include <memory>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/ScatteringFactorTable.hpp>

using namespace std;
using namespace diffpy::srreal;


class TestScatteringFactorTable : public CxxTest::TestSuite
{

    private:

        static const double mtol = 1.0e-4;
        ScatteringFactorTablePtr sftb;

    public:

        void test_factory()
        {
            ScatteringFactorTablePtr sfx0, sfx1, sfn0, sfn1;
            TS_ASSERT_THROWS(ScatteringFactorTable::createByType("invalid"),
                    invalid_argument);
            sfx0 = ScatteringFactorTable::createByType("periodictablexray");
            sfx1 = ScatteringFactorTable::createByType("X");
            TS_ASSERT(sfx0.get());
            TS_ASSERT(sfx1.get());
            TS_ASSERT_EQUALS(sfx0->type(), sfx1->type());
            TS_ASSERT(typeid(*sfx0) == typeid(*sfx1));
            sfn0 = ScatteringFactorTable::createByType(
                    "periodictableneutron");
            sfn1 = ScatteringFactorTable::createByType("N");
            TS_ASSERT(sfn0.get());
            TS_ASSERT(sfn1.get());
            TS_ASSERT_EQUALS(sfn0->type(), sfn1->type());
            TS_ASSERT(typeid(*sfn0) == typeid(*sfn1));
        }


        void test_setCustom()
        {
            sftb = ScatteringFactorTable::createByType("X");
            TS_ASSERT_DELTA(6.0, sftb->lookup("C"), 0.01);
            sftb->setCustom("C", 6.3);
            TS_ASSERT_THROWS(sftb->lookup("Ccustom"), invalid_argument);
            sftb->setCustom("Ccustom", 6.5);
            TS_ASSERT_EQUALS(6.5, sftb->lookup("Ccustom"));
            sftb->resetCustom("C");
            TS_ASSERT_EQUALS(6.5, sftb->lookup("Ccustom"));
            TS_ASSERT_DELTA(6.0, sftb->lookup("C"), 0.01);
            sftb->resetAll();
            TS_ASSERT_THROWS(sftb->lookup("Ccustom"), invalid_argument);
            TS_ASSERT_DELTA(6.0, sftb->lookup("C"), 0.01);
        }


        void test_periodictableXray()
        {
            sftb = ScatteringFactorTable::createByType("X");
            TS_ASSERT_DELTA(1.0, sftb->lookup("H"), 0.01);
            TS_ASSERT_DELTA(8.0, sftb->lookup("O"), 0.01);
            TS_ASSERT_DELTA(10.0, sftb->lookup("O2-"), 0.01);
            TS_ASSERT_DELTA(11.0, sftb->lookup("Na"), 0.01);
            TS_ASSERT_DELTA(10.0, sftb->lookup("Na+"), 0.01);
            TS_ASSERT_EQUALS(sftb->lookup("Na+"), sftb->lookup("Na1+"));
            TS_ASSERT_DELTA(74.0, sftb->lookup("W"), 0.04);
            TS_ASSERT_DELTA(88.0, sftb->lookup("Ra"), 0.04);
        }


        void test_periodictableNeutron()
        {
            sftb = ScatteringFactorTable::createByType("N");
            TS_ASSERT_DELTA(3.63, sftb->lookup("Na"), mtol);
            TS_ASSERT_DELTA(-3.37, sftb->lookup("Ti"), mtol);
            TS_ASSERT_DELTA(5.805, sftb->lookup("O"), mtol);
            TS_ASSERT_DELTA(6.6484, sftb->lookup("C"), mtol);
        }


        void test_ElectronNumber()
        {
            sftb = ScatteringFactorTable::createByType("electronnumber");
            TS_ASSERT_EQUALS(8.0, sftb->lookup("O"));
            TS_ASSERT_EQUALS(10.0, sftb->lookup("O2-"));
            TS_ASSERT_EQUALS(18.0, sftb->lookup("K+"));
            TS_ASSERT_EQUALS(18.0, sftb->lookup("K1+"));
            TS_ASSERT_EQUALS(68.0, sftb->lookup("W6+"));
            TS_ASSERT_THROWS(sftb->lookup("H4+"), invalid_argument);
            TS_ASSERT_THROWS(sftb->lookup("O0+"), invalid_argument);
        }

};

// End of file
