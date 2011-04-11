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
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/mathutils.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;


class TestScatteringFactorTable : public CxxTest::TestSuite
{

    private:

        static const double mtol = 1.0e-4;
        double meps;
        ScatteringFactorTablePtr msftb;

    public:

        void setUp()
        {
            meps = 10 * diffpy::mathutils::DOUBLE_EPS;
        }


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
            msftb = ScatteringFactorTable::createByType("X");
            TS_ASSERT_DELTA(6.0, msftb->lookup("C"), 0.01);
            msftb->setCustomFrom("C", "C", 6.3);
            TS_ASSERT_EQUALS(6.3, msftb->lookup("C"));
            msftb->setCustomFrom("C", "C", 6.4);
            TS_ASSERT_DELTA(6.4, msftb->lookup("C"), meps);
            TS_ASSERT_THROWS(msftb->lookup("Ccustom"), invalid_argument);
            msftb->setCustomFrom("Ccustom", "C", 6.5);
            TS_ASSERT_DELTA(6.5, msftb->lookup("Ccustom"), meps);
            msftb->resetCustom("C");
            TS_ASSERT_DELTA(6.5, msftb->lookup("Ccustom"), meps);
            TS_ASSERT_DELTA(6.0, msftb->lookup("C"), 0.01);
            msftb->resetAll();
            TS_ASSERT_THROWS(msftb->lookup("Ccustom"), invalid_argument);
            TS_ASSERT_DELTA(6.0, msftb->lookup("C"), 0.01);
        }


        void test_getCustomSymbols()
        {
            msftb = ScatteringFactorTable::createByType("X");
            TS_ASSERT(msftb->getCustomSymbols().empty());
            msftb->setCustomFrom("C", "C", 6.1);
            TS_ASSERT_EQUALS(1u, msftb->getCustomSymbols().size());
            TS_ASSERT_EQUALS(1u, msftb->getCustomSymbols().count("C"));
            msftb->setCustomFrom("C", "C", 6.3);
            TS_ASSERT_EQUALS(1u, msftb->getCustomSymbols().size());
            TS_ASSERT_EQUALS(1u, msftb->getCustomSymbols().count("C"));
            TS_ASSERT_EQUALS(6.3, msftb->lookup("C"));
            ScatteringFactorTablePtr sftb1 = msftb->clone();
            msftb->resetCustom("C");
            TS_ASSERT(msftb->getCustomSymbols().empty());
            TS_ASSERT_EQUALS(1u, sftb1->getCustomSymbols().size());
            TS_ASSERT_EQUALS(1u, sftb1->getCustomSymbols().count("C"));
            TS_ASSERT_EQUALS(6.3, sftb1->lookup("C"));
            sftb1->resetAll();
            TS_ASSERT(msftb->getCustomSymbols().empty());
        }


        void test_periodictableXray()
        {
            msftb = ScatteringFactorTable::createByType("X");
            TS_ASSERT_DELTA(1.0, msftb->lookup("H"), 0.01);
            TS_ASSERT_DELTA(8.0, msftb->lookup("O"), 0.01);
            TS_ASSERT_DELTA(10.0, msftb->lookup("O2-"), 0.01);
            TS_ASSERT_DELTA(11.0, msftb->lookup("Na"), 0.01);
            TS_ASSERT_DELTA(10.0, msftb->lookup("Na+"), 0.01);
            TS_ASSERT_EQUALS(msftb->lookup("Na+"), msftb->lookup("Na1+"));
            TS_ASSERT_DELTA(74.0, msftb->lookup("W"), 0.04);
            TS_ASSERT_DELTA(88.0, msftb->lookup("Ra"), 0.04);
        }


        void test_periodictableElectron()
        {
            using diffpy::mathutils::DOUBLE_MAX;
            msftb = ScatteringFactorTable::createByType("E");
            TS_ASSERT_EQUALS(DOUBLE_MAX, msftb->lookup("H"));
            TS_ASSERT_EQUALS(DOUBLE_MAX, msftb->lookup("Ra"));
            TS_ASSERT_DELTA(3.42104, msftb->lookup("Na", 1), 1e-5);
            TS_ASSERT_DELTA(1.34868, msftb->lookup("Na", 3), 1e-5);
            TS_ASSERT_DELTA(0.832158, msftb->lookup("Na", 5), 1e-5);
            TS_ASSERT_THROWS(msftb->lookup("H4+"), invalid_argument);
            TS_ASSERT_THROWS(msftb->lookup("H4+", 3), invalid_argument);
        }


        void test_periodictableNeutron()
        {
            msftb = ScatteringFactorTable::createByType("N");
            TS_ASSERT_DELTA(3.63, msftb->lookup("Na"), mtol);
            TS_ASSERT_DELTA(-3.37, msftb->lookup("Ti"), mtol);
            TS_ASSERT_DELTA(5.805, msftb->lookup("O"), mtol);
            TS_ASSERT_DELTA(6.6484, msftb->lookup("C"), mtol);
        }


        void test_ElectronNumber()
        {
            msftb = ScatteringFactorTable::createByType("electronnumber");
            TS_ASSERT_EQUALS(8.0, msftb->lookup("O"));
            TS_ASSERT_EQUALS(10.0, msftb->lookup("O2-"));
            TS_ASSERT_EQUALS(18.0, msftb->lookup("K+"));
            TS_ASSERT_EQUALS(18.0, msftb->lookup("K1+"));
            TS_ASSERT_EQUALS(68.0, msftb->lookup("W6+"));
            TS_ASSERT_THROWS(msftb->lookup("H4+"), invalid_argument);
            TS_ASSERT_THROWS(msftb->lookup("O0+"), invalid_argument);
        }


        void test_serialization()
        {
            ScatteringFactorTablePtr sftb1;
            msftb = ScatteringFactorTable::createByType("electronnumber");
            msftb->setCustomFrom("H", "H", 1.23);
            sftb1 = dumpandload(msftb);
            TS_ASSERT_EQUALS(string("electronnumber"), sftb1->type());
            TS_ASSERT_EQUALS(1.23, sftb1->lookup("H"));
            TS_ASSERT_EQUALS(1u, sftb1->getCustomSymbols().size());
            sftb1 = dumpandload(ScatteringFactorTable::createByType("N"));
            TS_ASSERT_EQUALS(string("periodictableneutron"), sftb1->type());
            sftb1 = dumpandload(ScatteringFactorTable::createByType("X"));
            TS_ASSERT_EQUALS(string("periodictablexray"), sftb1->type());
        }

};

// End of file
