/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class TestPDFCalculator -- unit tests for PDFCalculator class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/QResolutionEnvelope.hpp>
#include <diffpy/serialization.hpp>
#include "test_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

class TestPDFCalculator : public CxxTest::TestSuite
{
    private:

        boost::shared_ptr<PDFCalculator> mpdfc;
        StructureAdapterPtr memptystru;
        double meps;
        double mepsdb;

    public:

        void setUp()
        {
            memptystru = emptyStructureAdapter();
            meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
            mepsdb = 10 * diffpy::mathutils::DOUBLE_EPS;
            mpdfc.reset(new PDFCalculator);
        }


        void test_setPeakWidthModel()
        {
            const PeakWidthModelPtr& jpw0 = mpdfc->getPeakWidthModel();
            TS_ASSERT_EQUALS(0.0, jpw0->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(0.0, jpw0->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(0.0, jpw0->getDoubleAttr("qbroad"));
            TS_ASSERT_EQUALS(0.0, jpw0->getDoubleAttr("qbroad_new"));
            JeongPeakWidth jpw;
            jpw.setDelta1(1.0);
            jpw.setDelta2(2.0);
            jpw.setQbroad(3.0);
            jpw.setQbroad_new(4.0);
            mpdfc->setPeakWidthModel(jpw.clone());
            const PeakWidthModelPtr& jpw1 = mpdfc->getPeakWidthModel();
            TS_ASSERT_EQUALS(1.0, jpw1->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(2.0, jpw1->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(3.0, jpw1->getDoubleAttr("qbroad"));
            TS_ASSERT_EQUALS(4.0, jpw1->getDoubleAttr("qbroad_new"));
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
            const double qstep0 = 100 * M_PI / 1024;
            const double qstep1 = 100 * M_PI / 2048;
            const double qstep2 = 100 * M_PI / 4096;
            const double qstep3 = M_PI / (1024 * 0.03);
            TS_ASSERT_DELTA(qstep0, mpdfc->getQstep(), meps);
            mpdfc->setRmax(20);
            TS_ASSERT_DELTA(qstep1, mpdfc->getQstep(), meps);
            mpdfc->setQmax(10);
            TS_ASSERT_DELTA(qstep2, mpdfc->getQstep(), meps);
            mpdfc->setRstep(0.03);
            TS_ASSERT_DELTA(qstep3, mpdfc->getQstep(), meps);
            // test qstep after evaluation of some structure.
            StructureAdapterPtr catio3;
            catio3 = loadTestPeriodicStructure("CaTiO3.stru");
            mpdfc->eval(catio3);
            TS_ASSERT_DELTA(qstep3, mpdfc->getQstep(), meps);
            // now qstep is not updated because of non-trivial structure.
            mpdfc->setRstep(0.01);
            TS_ASSERT_DIFFERS(qstep2, mpdfc->getQstep());
            // qstep is correct after calculation.
            mpdfc->eval();
            TS_ASSERT_DELTA(qstep2, mpdfc->getQstep(), meps);
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
            TS_ASSERT_DELTA(5.0, mpdfc->getExtendedRmin(), mepsdb);
        }


        void test_getExtendedRmax()
        {
            // empty structure should not extend the grid at all.
            TS_ASSERT_DELTA(10.0, mpdfc->getExtendedRmax(), mepsdb);
            mpdfc->setRmax(7);
            mpdfc->eval(memptystru);
            TS_ASSERT_DELTA(7.0, mpdfc->getExtendedRmax(), mepsdb);
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
