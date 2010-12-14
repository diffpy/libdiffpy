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
* class TestStepCutEnvelope -- unit tests for StepCutEnvelope class
*
* $Id$
*
*****************************************************************************/

#include <memory>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/PDFEnvelope.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

class TestStepCutEnvelope : public CxxTest::TestSuite
{
    private:

        PDFEnvelopePtr menvelope;

    public:

        void setUp()
        {
            menvelope = PDFEnvelope::createByType("stepcut");
        }


        void test_create()
        {
            TS_ASSERT_EQUALS(0.0, menvelope->getDoubleAttr("stepcut"));
            menvelope->setDoubleAttr("stepcut", 13.0);
            TS_ASSERT_EQUALS(13.0, menvelope->getDoubleAttr("stepcut"));
            PDFEnvelopePtr e1(menvelope->create());
            TS_ASSERT_EQUALS(0.0, e1->getDoubleAttr("stepcut"));
        }


        void test_copy()
        {
            menvelope->setDoubleAttr("stepcut", 13.0);
            TS_ASSERT_EQUALS(13.0, menvelope->getDoubleAttr("stepcut"));
            PDFEnvelopePtr e1(menvelope->clone());
            TS_ASSERT_EQUALS(13.0, e1->getDoubleAttr("stepcut"));
        }


        void test_type()
        {
            TS_ASSERT_EQUALS("stepcut", menvelope->type());
        }


        void test_parentheses_operator()
        {
            const PDFEnvelope& fne = *menvelope;
            TS_ASSERT_EQUALS(1.0, fne(-1.0));
            TS_ASSERT_EQUALS(1.0, fne(0.0));
            TS_ASSERT_EQUALS(1.0, fne(+1.0));
            menvelope->setDoubleAttr("stepcut", 1.0);
            TS_ASSERT_EQUALS(1.0, fne(-1.0));
            TS_ASSERT_EQUALS(1.0, fne(0.0));
            TS_ASSERT_EQUALS(1.0, fne(+1.0));
            TS_ASSERT_EQUALS(0.0, fne(+1.0001));
            TS_ASSERT_EQUALS(0.0, fne(+2.0));
        }


        void test_serialization()
        {
            menvelope->setDoubleAttr("stepcut", 13.1);
            PDFEnvelopePtr e1 = dumpandload(menvelope);
            TS_ASSERT_EQUALS(string("stepcut"), e1->type());
            TS_ASSERT_EQUALS(13.1, e1->getDoubleAttr("stepcut"));
        }

};  // class TestDiffPyStructureBondGenerator

// End of file
