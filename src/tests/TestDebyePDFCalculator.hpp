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
* class TestDebyePDFCalculator -- unit tests for DebyePDFCalculator class
*
* $Id$
*
*****************************************************************************/

#include <memory>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/QResolutionEnvelope.hpp>
#include <diffpy/srreal/VR3Structure.hpp>
#include <diffpy/serialization.hpp>

using namespace std;
using namespace diffpy::srreal;

class TestDebyePDFCalculator : public CxxTest::TestSuite
{
    private:

        boost::shared_ptr<DebyePDFCalculator> mpdfc;
        VR3Structure memptystru;
        double meps;

    public:

        void setUp()
        {
            meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
            mpdfc.reset(new DebyePDFCalculator);
        }


        void test_access_Envelopes()
        {
            TS_ASSERT_EQUALS(2u, mpdfc->usedEnvelopeTypes().size());
            TS_ASSERT_EQUALS(1.0, mpdfc->getDoubleAttr("scale"));
            TS_ASSERT_EQUALS(0.0, mpdfc->getDoubleAttr("qdamp"));
            mpdfc->setDoubleAttr("scale", 3.0);
            TS_ASSERT_EQUALS(3.0, mpdfc->getDoubleAttr("scale"));
            mpdfc->addEnvelopeByType("scale");
            TS_ASSERT_EQUALS(1.0, mpdfc->getDoubleAttr("scale"));
            QResolutionEnvelope qdamp4;
            qdamp4.setQdamp(4);
            mpdfc->addEnvelope(qdamp4.clone());
            TS_ASSERT_EQUALS(4.0, mpdfc->getDoubleAttr("qdamp"));
            TS_ASSERT_THROWS(mpdfc->addEnvelopeByType("invalid"), logic_error);
        }


        void test_getPDF()
        {
            QuantityType pdf;
            TS_ASSERT_EQUALS(1000u, mpdfc->getPDF().size());
            mpdfc->setRmin(2.0);
            mpdfc->setRmax(0.0);
            mpdfc->eval(memptystru);
            pdf = mpdfc->getPDF();
            TS_ASSERT(mpdfc->getPDF().empty());
            mpdfc->setRmax(2.0);
            mpdfc->eval(memptystru);
            pdf = mpdfc->getPDF();
            TS_ASSERT(mpdfc->getPDF().empty());
            mpdfc->setRmax(2.01001);
            mpdfc->eval(memptystru);
            pdf = mpdfc->getPDF();
            TS_ASSERT_EQUALS(2u, pdf.size());
        }


        void test_getRDF()
        {
            // getRDF is not yet implemented and just returns empty array.
            QuantityType rdf = mpdfc->getRDF();
            TS_ASSERT(rdf.empty());

        }


        void test_getF()
        {
            QuantityType fq = mpdfc->getF();
            TS_ASSERT_EQUALS(92u, fq.size());
            TS_ASSERT_EQUALS(0.0, *min_element(fq.begin(), fq.end()));
            TS_ASSERT_EQUALS(0.0, *max_element(fq.begin(), fq.end()));
        }


        void test_getQgrid()
        {
            TS_ASSERT_EQUALS(92u, mpdfc->getQgrid().size());
        }


        void test_getQmin()
        {
            TS_ASSERT_EQUALS(0.0, mpdfc->getQmin());
        }


        void test_getQmax()
        {
            TS_ASSERT_DELTA(25.0, mpdfc->getQmax(), meps);
        }


        void test_getQstep()
        {
            const double rcalchi = ceil((10.0 + 12*M_PI/25.0) / 0.01) * 0.01;
            const double qstep = M_PI / rcalchi;
            TS_ASSERT_DELTA(qstep, mpdfc->getQstep(), meps);
        }


        void test_getRgrid()
        {
            QuantityType rgrid0 = mpdfc->getRgrid();
            TS_ASSERT_EQUALS(rgrid0, mpdfc->getRgrid());
            mpdfc->setRmin(5);
            mpdfc->setRmax(4);
            TS_ASSERT(mpdfc->getRgrid().empty());
        }


        void test_serialization()
        {
            // build customized DebyePDFCalculator
            mpdfc->setPeakWidthModelByType("constant");
            mpdfc->setDoubleAttr("width", 0.123);
            mpdfc->setDoubleAttr("debyeprecision", 0.00011);
            mpdfc->setScatteringFactorTableByType("electronnumber");
            mpdfc->getScatteringFactorTable()->setCustom("H", 1.1);
            // dump it to string
            stringstream storage(ios::in | ios::out | ios::binary);
            diffpy::serialization::oarchive oa(storage, ios::binary);
            oa << mpdfc;
            diffpy::serialization::iarchive ia(storage, ios::binary);
            boost::shared_ptr<DebyePDFCalculator> pdfc1;
            ia >> pdfc1;
            TS_ASSERT_DIFFERS(pdfc1.get(), mpdfc.get());
            TS_ASSERT_EQUALS(string("constant"),
                    pdfc1->getPeakWidthModel()->type());
            TS_ASSERT_EQUALS(0.123, pdfc1->getDoubleAttr("width"));
            TS_ASSERT_EQUALS(0.00011, pdfc1->getDoubleAttr("debyeprecision"));
            TS_ASSERT_EQUALS(string("electronnumber"),
                    pdfc1->getScatteringFactorTable()->type());
            TS_ASSERT_EQUALS(1.1,
                    pdfc1->getScatteringFactorTable()->lookup("H"));
        }

};  // class TestDebyePDFCalculator

// End of file
