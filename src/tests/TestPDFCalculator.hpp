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
* class TestPDFCalculator -- unit tests for PDFCalculator class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/VR3Structure.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/QResolutionEnvelope.hpp>
#include <diffpy/serialization.hpp>

using namespace std;
using namespace diffpy::srreal;

class TestPDFCalculator : public CxxTest::TestSuite
{
    private:

        boost::shared_ptr<PDFCalculator> mpdfc;
        VR3Structure memptystru;
        double meps;

    public:

        void setUp()
        {
            meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
            mpdfc.reset(new PDFCalculator);
        }


        void test_setPeakWidthModel()
        {
            const PeakWidthModelPtr& jpw0 = mpdfc->getPeakWidthModel();
            TS_ASSERT_EQUALS(0.0, jpw0->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(0.0, jpw0->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(0.0, jpw0->getDoubleAttr("qbroad"));
            JeongPeakWidth jpw;
            jpw.setDelta1(1.0);
            jpw.setDelta2(2.0);
            jpw.setQbroad(3.0);
            mpdfc->setPeakWidthModel(jpw.clone());
            const PeakWidthModelPtr& jpw1 = mpdfc->getPeakWidthModel();
            TS_ASSERT_EQUALS(1.0, jpw1->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(2.0, jpw1->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(3.0, jpw1->getDoubleAttr("qbroad"));
        }


        void test_getPeakWidthModel()
        {
            string tp = "jeong";
            TS_ASSERT_EQUALS(tp, mpdfc->getPeakWidthModel()->type());
            PeakWidthModelPtr pwm(PeakWidthModel::createByType("debye-waller"));
            mpdfc->setPeakWidthModel(pwm);
            tp = "debye-waller";
            TS_ASSERT_EQUALS(tp, mpdfc->getPeakWidthModel()->type());
            mpdfc->setPeakWidthModel(ConstantPeakWidth().clone());
            tp = "constant";
            TS_ASSERT_EQUALS(tp, mpdfc->getPeakWidthModel()->type());
        }


        void test_access_PeakProfile()
        {
            PeakProfilePtr pkf(mpdfc->getPeakProfile()->clone());
            pkf->setPrecision(1.1);
            mpdfc->setPeakProfile(pkf);
            TS_ASSERT_EQUALS(1.1, mpdfc->getPeakProfile()->getPrecision());
            TS_ASSERT_EQUALS(1.1, mpdfc->getDoubleAttr("peakprecision"));
            mpdfc->setDoubleAttr("peakprecision", 0.2);
            TS_ASSERT_EQUALS(0.2, mpdfc->getDoubleAttr("peakprecision"));
            TS_ASSERT_EQUALS(pkf->type(), mpdfc->getPeakProfile()->type());
            mpdfc->setPeakProfileByType("gaussian");
            TS_ASSERT_EQUALS(0.2, mpdfc->getDoubleAttr("peakprecision"));
            TS_ASSERT_THROWS(mpdfc->setPeakProfileByType("invalid"),
                    logic_error);
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
            pdf = mpdfc->getPDF();
            TS_ASSERT_EQUALS(1000u, pdf.size());
            mpdfc->setRmin(2.0);
            mpdfc->setRmax(0.0);
            mpdfc->eval(memptystru);
            pdf = mpdfc->getPDF();
            TS_ASSERT(pdf.empty());
            mpdfc->setRmax(2.0);
            mpdfc->eval(memptystru);
            pdf = mpdfc->getPDF();
            TS_ASSERT(pdf.empty());
            mpdfc->setRmax(2.01001);
            mpdfc->eval(memptystru);
            pdf = mpdfc->getPDF();
            TS_ASSERT_EQUALS(2u, pdf.size());
            TS_ASSERT_EQUALS(0.0, pdf[0]);
            TS_ASSERT_EQUALS(0.0, pdf[1]);
        }


        void test_getRDF()
        {
            QuantityType rdf = mpdfc->getRDF();
            TS_ASSERT_EQUALS(1000u, rdf.size());
            TS_ASSERT_EQUALS(0.0, *min_element(rdf.begin(), rdf.end()));
            TS_ASSERT_EQUALS(0.0, *max_element(rdf.begin(), rdf.end()));
            mpdfc->eval(memptystru);
            rdf = mpdfc->getRDF();
            TS_ASSERT_EQUALS(1000u, rdf.size());
            TS_ASSERT_EQUALS(0.0, *min_element(rdf.begin(), rdf.end()));
            TS_ASSERT_EQUALS(0.0, *max_element(rdf.begin(), rdf.end()));
            mpdfc->setRmin(2.0);
            mpdfc->setRmax(0.0);
            mpdfc->eval(memptystru);
            rdf = mpdfc->getRDF();
            TS_ASSERT(rdf.empty());
            mpdfc->setRmax(2.01001);
            mpdfc->eval(memptystru);
            rdf = mpdfc->getRDF();
            TS_ASSERT_EQUALS(2u, rdf.size());
            TS_ASSERT_EQUALS(0.0, rdf[0]);
            TS_ASSERT_EQUALS(0.0, rdf[1]);

        }


        void test_getF()
        {
            QuantityType fq = mpdfc->getF();
            TS_ASSERT_EQUALS(1024u, fq.size());
            TS_ASSERT_EQUALS(0.0, *min_element(fq.begin(), fq.end()));
            TS_ASSERT_EQUALS(0.0, *max_element(fq.begin(), fq.end()));
        }


        void test_getQgrid()
        {
            TS_ASSERT_EQUALS(1024u, mpdfc->getQgrid().size());
        }


        void test_getQmin()
        {
            TS_ASSERT_EQUALS(0.0, mpdfc->getQmin());
        }


        void test_getQmax()
        {
            TS_ASSERT_DELTA(100 * M_PI, mpdfc->getQmax(), meps);
        }


        void test_getQstep()
        {
            TS_ASSERT_DELTA(100 * M_PI / 1024, mpdfc->getQstep(), meps);
        }


        void test_getRgrid()
        {
            QuantityType rgrid0 = mpdfc->getRgrid();
            TS_ASSERT_EQUALS(rgrid0, mpdfc->getRgrid());
            mpdfc->setRmin(5);
            mpdfc->setRmax(4);
            TS_ASSERT(mpdfc->getRgrid().empty());
        }


        void test_getExtendedRmin()
        {
            // empty structure should not extend the grid at all.
            TS_ASSERT_EQUALS(0.0, mpdfc->getExtendedRmin());
            mpdfc->setRmin(5);
            mpdfc->eval(memptystru);
            TS_ASSERT_EQUALS(5.0, mpdfc->getExtendedRmin());
        }


        void test_getExtendedRmax()
        {
            // empty structure should not extend the grid at all.
            TS_ASSERT_EQUALS(10.0, mpdfc->getExtendedRmax());
            mpdfc->setRmax(7);
            mpdfc->eval(memptystru);
            TS_ASSERT_EQUALS(7.0, mpdfc->getExtendedRmax());
        }


        void test_serialization()
        {
            // build customized PDFCalculator
            mpdfc->setPeakWidthModelByType("constant");
            mpdfc->setDoubleAttr("width", 0.123);
            mpdfc->setPeakProfileByType("gaussian");
            mpdfc->setDoubleAttr("peakprecision", 0.011);
            mpdfc->setScatteringFactorTableByType("electronnumber");
            mpdfc->getScatteringFactorTable()->setCustomAs("H", "H", 1.1);
            // dump it to string
            stringstream storage(ios::in | ios::out | ios::binary);
            diffpy::serialization::oarchive oa(storage, ios::binary);
            oa << mpdfc;
            diffpy::serialization::iarchive ia(storage, ios::binary);
            boost::shared_ptr<PDFCalculator> pdfc1;
            ia >> pdfc1;
            TS_ASSERT_DIFFERS(pdfc1.get(), mpdfc.get());
            TS_ASSERT_EQUALS(string("constant"),
                    pdfc1->getPeakWidthModel()->type());
            TS_ASSERT_EQUALS(0.123, pdfc1->getDoubleAttr("width"));
            TS_ASSERT_EQUALS(0.011, pdfc1->getDoubleAttr("peakprecision"));
            TS_ASSERT_EQUALS(string("electronnumber"),
                    pdfc1->getScatteringFactorTable()->type());
            TS_ASSERT_EQUALS(1.1,
                    pdfc1->getScatteringFactorTable()->lookup("H"));
        }

};  // class TestPDFCalculator

// End of file
