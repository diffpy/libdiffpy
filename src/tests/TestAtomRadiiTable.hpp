/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestAtomRadiiTable -- unit tests for the AtomRadiiTable class
*
*****************************************************************************/

#include <stdexcept>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/AtomRadiiTable.hpp>
#include <diffpy/srreal/ConstantRadiiTable.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;


class TestAtomRadiiTable : public CxxTest::TestSuite
{

    private:

        AtomRadiiTablePtr mrtb;

    public:

        void setUp()
        {
            mrtb = AtomRadiiTable::createByType("constant");
        }


        void test_lookup()
        {
            TS_ASSERT_EQUALS(0.0, mrtb->lookup("everything"));
        }


        void test_setCustom()
        {
            mrtb->setCustom("C", 1.23);
            TS_ASSERT_EQUALS(1.23, mrtb->lookup("C"));
            mrtb->setCustom("C", 1.5);
            TS_ASSERT_EQUALS(1.5, mrtb->lookup("C"));
            mrtb->resetAll();
            TS_ASSERT_EQUALS(0.0, mrtb->lookup("C"));
        }


        void test_fromString()
        {
            TS_ASSERT_EQUALS(0u, mrtb->getAllCustom().size());
            mrtb->fromString("A:1.1, B:1.2,    C:  1.3  \t");
            TS_ASSERT_EQUALS(3u, mrtb->getAllCustom().size());
            TS_ASSERT_EQUALS(1.1, mrtb->lookup("A"));
            TS_ASSERT_EQUALS(1.2, mrtb->lookup("B"));
            TS_ASSERT_EQUALS(1.3, mrtb->lookup("C"));
            TS_ASSERT_THROWS(mrtb->fromString("X12"), invalid_argument);
            TS_ASSERT_EQUALS(3u, mrtb->getAllCustom().size());
            mrtb->fromString("D:4.3");
            TS_ASSERT_EQUALS(4.3, mrtb->lookup("D"));
            TS_ASSERT_EQUALS(4u, mrtb->getAllCustom().size());
        }


        void test_resetCustom()
        {
            mrtb->fromString("A:1.1, B:1.2,    C:  1.3");
            TS_ASSERT_EQUALS(3u, mrtb->getAllCustom().size());
            mrtb->resetCustom("B");
            TS_ASSERT_EQUALS(1.1, mrtb->lookup("A"));
            TS_ASSERT_EQUALS(1.3, mrtb->lookup("C"));
        }


        void test_resetAll()
        {
            mrtb->fromString("A:1.1, B:1.2,    C:  1.3");
            TS_ASSERT_EQUALS(3u, mrtb->getAllCustom().size());
            mrtb->resetAll();
            TS_ASSERT(mrtb->getAllCustom().empty());
        }


        void test_toString()
        {
            TS_ASSERT(mrtb->toString().empty());
            mrtb->fromString("A:1.1,    C:  1.3, B:1.2");
            TS_ASSERT_EQUALS(string("A:1.1,B:1.2,C:1.3"), mrtb->toString());
            TS_ASSERT_EQUALS(string("A:1.1, B:1.2, C:1.3"), mrtb->toString(", "));
            mrtb->resetAll();
            TS_ASSERT(mrtb->toString().empty());
        }


        void test_setDefault()
        {
            ConstantRadiiTable* crtb =
                dynamic_cast<ConstantRadiiTable*>(mrtb.get());
            crtb->setDefault(7.1);
            TS_ASSERT_EQUALS(7.1, mrtb->lookup("everything"));
            crtb->setDefault(3.1);
            TS_ASSERT_EQUALS(3.1, mrtb->lookup("Na"));
        }


        void test_serialization()
        {
            AtomRadiiTablePtr rtb1;
            mrtb->setCustom("H", 1.23);
            ConstantRadiiTable* crtb =
                dynamic_cast<ConstantRadiiTable*>(mrtb.get());
            crtb->setDefault(0.5);
            rtb1 = dumpandload(mrtb);
            TS_ASSERT_EQUALS(1.23, rtb1->lookup("H"));
            TS_ASSERT_EQUALS(1u, rtb1->getAllCustom().size());
            TS_ASSERT_DIFFERS(mrtb.get(), rtb1.get());
            ConstantRadiiTable* crtb1 =
                dynamic_cast<ConstantRadiiTable*>(rtb1.get());
            TS_ASSERT_EQUALS(0.5, crtb1->getDefault());
        }

};

// End of file
