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
* class TestSphericalShapeEnvelope -- unit tests for SphericalShapeEnvelope
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/PDFEnvelope.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

class TestSphericalShapeEnvelope : public CxxTest::TestSuite
{
    private:

        PDFEnvelopePtr menvelope;

    public:

        void setUp()
        {
            menvelope = PDFEnvelope::createByType("sphericalshape");
        }


        void test_create()
        {
            TS_ASSERT_EQUALS(0.0, menvelope->getDoubleAttr("spdiameter"));
            menvelope->setDoubleAttr("spdiameter", 13.0);
            TS_ASSERT_EQUALS(13.0, menvelope->getDoubleAttr("spdiameter"));
            PDFEnvelopePtr e1 = menvelope->create();
            TS_ASSERT_EQUALS(0.0, e1->getDoubleAttr("spdiameter"));
        }


        void test_copy()
        {
            menvelope->setDoubleAttr("spdiameter", 13.0);
            TS_ASSERT_EQUALS(13.0, menvelope->getDoubleAttr("spdiameter"));
            PDFEnvelopePtr e1 = menvelope->clone();
            TS_ASSERT_EQUALS(13.0, e1->getDoubleAttr("spdiameter"));
        }


        void test_type()
        {
            TS_ASSERT_EQUALS("sphericalshape", menvelope->type());
        }


        void test_parentheses_operator()
        {
            const PDFEnvelope& fne = *menvelope;
            TS_ASSERT_EQUALS(1.0, fne(0.0));
            TS_ASSERT_EQUALS(1.0, fne(100.0));
            menvelope->setDoubleAttr("spdiameter", 10.0);
            TS_ASSERT_EQUALS(1.0, fne(0.0));
            TS_ASSERT_EQUALS(0.0, fne(10.0));
            TS_ASSERT_EQUALS(0.0, fne(100));
            TS_ASSERT(fne(9.99) > 0.0);
            TS_ASSERT_DELTA(0.3125, fne(5), 1e-8);
        }


        void test_serialization()
        {
            menvelope->setDoubleAttr("spdiameter", 13.1);
            PDFEnvelopePtr e1 = dumpandload(menvelope);
            TS_ASSERT_EQUALS(string("sphericalshape"), e1->type());
            TS_ASSERT_EQUALS(13.1, e1->getDoubleAttr("spdiameter"));
        }

};  // class TestSphericalShapeEnvelope

// End of file
