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
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <boost/make_shared.hpp>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/QResolutionEnvelope.hpp>
#include <diffpy/serialization.hpp>

using namespace std;
using namespace diffpy::srreal;

class TestDebyePDFCalculator : public CxxTest::TestSuite
{
    private:

        boost::shared_ptr<DebyePDFCalculator> mpdfc;
        AtomicStructureAdapterPtr memptystru;
        AtomicStructureAdapterPtr mstru10;
        AtomicStructureAdapterPtr mstru10d1;
        AtomicStructureAdapterPtr mstru10r;
        AtomicStructureAdapterPtr mstru9;
        diffpy::mathutils::EpsilonEqual allclose;
        double meps;

    public:

        void setUp()
        {
            const int SZ = 10;
            meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
            mpdfc.reset(new DebyePDFCalculator);
            memptystru = boost::make_shared<AtomicStructureAdapter>();
            mstru10 = boost::make_shared<AtomicStructureAdapter>();
            Atom ai;
            ai.atomtype = "C";
            ai.cartesianuij = R3::identity();
            ai.cartesianuij(0, 0) = ai.cartesianuij(1, 1) =
                ai.cartesianuij(2, 2) = 0.004;
            for (int i = 0; i < SZ; ++i)
            {
                ai.cartesianposition[0] = i;
                mstru10->append(ai);
            }
            mstru10d1 = boost::make_shared<AtomicStructureAdapter>(*mstru10);
            (*mstru10d1)[0].atomtype = "Au";
            mstru10r = boost::make_shared<AtomicStructureAdapter>();
            mstru10r->assign(mstru10->rbegin(), mstru10->rend());
            mstru9 = boost::make_shared<AtomicStructureAdapter>(*mstru10);
            mstru9->remove(9);
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
            TS_ASSERT_EQUALS(1000u, rdf.size());
            TS_ASSERT_EQUALS(0.0, *std::min_element(rdf.begin(), rdf.end()));
            TS_ASSERT_EQUALS(0.0, *std::max_element(rdf.begin(), rdf.end()));
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
            mpdfc->setQmin(1.0);
            const DebyePDFCalculator& dpc = *mpdfc;
            TS_ASSERT_EQUALS(1.0, dpc.getQmin());
            TS_ASSERT_EQUALS(0.0, dpc.BaseDebyeSum::getQmin());
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
            mpdfc->getScatteringFactorTable()->setCustomAs("H", "H", 1.1);
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


        void test_DBPDF_change_atom()
        {
            mpdfc->setQmin(1.0);
            DebyePDFCalculator pdfcb = *mpdfc;
            DebyePDFCalculator pdfco = *mpdfc;
            TS_ASSERT_EQUALS(1.0, pdfcb.getQmin());
            pdfcb.setEvaluatorType(BASIC);
            pdfco.setEvaluatorType(OPTIMIZED);
            TS_ASSERT_EQUALS(NONE, pdfcb.getEvaluatorTypeUsed());
            TS_ASSERT_EQUALS(NONE, pdfco.getEvaluatorTypeUsed());
            pdfcb.eval(mstru10);
            pdfco.eval(mstru10);
            QuantityType gb = pdfcb.getPDF();
            QuantityType go = pdfco.getPDF();
            TS_ASSERT(!gb.empty());
            TS_ASSERT_EQUALS(gb, go);
            int cnonzero = count_if(gb.begin(), gb.end(),
                    bind1st(not_equal_to<double>(), 0.0));
            TS_ASSERT(cnonzero);
            TS_ASSERT_EQUALS(BASIC, pdfcb.getEvaluatorType());
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.getEvaluatorType());
            // first call of pdfco should use the BASIC evaluation
            TS_ASSERT_EQUALS(BASIC, pdfcb.getEvaluatorTypeUsed());
            TS_ASSERT_EQUALS(BASIC, pdfco.getEvaluatorTypeUsed());
            // test second call on the same structure
            pdfco.eval(mstru10);
            go = pdfco.getPDF();
            TS_ASSERT_EQUALS(gb, go);
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.getEvaluatorTypeUsed());
            // test structure with one different atom
            pdfcb.eval(mstru10d1);
            pdfco.eval(mstru10d1);
            QuantityType gb1 = pdfcb.getPDF();
            QuantityType go1 = pdfco.getPDF();
            TS_ASSERT(!allclose(gb, gb1));
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.getEvaluatorTypeUsed());
            TS_ASSERT(allclose(gb1, go1));
            // change position of 1 atom
            AtomicStructureAdapterPtr stru10d1s =
                boost::make_shared<AtomicStructureAdapter>(*mstru10d1);
            (*stru10d1s)[0].cartesianposition[1] = 0.5;
            pdfcb.eval(stru10d1s);
            pdfco.eval(stru10d1s);
            QuantityType gb2 = pdfcb.getPDF();
            QuantityType go2 = pdfco.getPDF();
            TS_ASSERT(!allclose(gb1, gb2));
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.getEvaluatorTypeUsed());
            TS_ASSERT(allclose(gb2, go2));
        }


        void test_DBPDF_reverse_atoms()
        {
            mpdfc->setQmin(1.1);
            DebyePDFCalculator pdfcb = *mpdfc;
            DebyePDFCalculator pdfco = *mpdfc;
            pdfcb.setEvaluatorType(BASIC);
            pdfco.setEvaluatorType(OPTIMIZED);
            pdfcb.eval(mstru10);
            pdfco.eval(mstru10);
            pdfco.eval(mstru10r);
            QuantityType gb = pdfcb.getPDF();
            QuantityType go = pdfco.getPDF();
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.getEvaluatorTypeUsed());
            TS_ASSERT(allclose(gb, go));
        }


        void test_DBPDF_remove_atom()
        {
            mpdfc->setQmin(1.2);
            DebyePDFCalculator pdfcb = *mpdfc;
            DebyePDFCalculator pdfco = *mpdfc;
            pdfcb.setEvaluatorType(BASIC);
            pdfco.setEvaluatorType(OPTIMIZED);
            pdfcb.eval(mstru9);
            pdfco.eval(mstru10);
            pdfco.eval(mstru9);
            QuantityType gb = pdfcb.getPDF();
            QuantityType go = pdfco.getPDF();
            TS_ASSERT_EQUALS(OPTIMIZED, pdfco.getEvaluatorTypeUsed());
            TS_ASSERT(allclose(gb, go));
        }


        void test_DBPDF_qmin_click()
        {
            mpdfc->setEvaluatorType(OPTIMIZED);
            TS_ASSERT_EQUALS(NONE, mpdfc->getEvaluatorTypeUsed());
            mpdfc->eval(mstru10);
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfc->getEvaluatorTypeUsed());
            mpdfc->setQmin(4);
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfc->getEvaluatorTypeUsed());
            mpdfc->setQmin(0);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfc->getEvaluatorTypeUsed());
        }


        void test_DBPDF_qmax_click()
        {
            mpdfc->setEvaluatorType(OPTIMIZED);
            mpdfc->eval(mstru10);
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfc->getEvaluatorTypeUsed());
            mpdfc->setQmax(mpdfc->getQmax());
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfc->getEvaluatorTypeUsed());
            mpdfc->setQmax(mpdfc->getQmax() - 1.0);
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(BASIC, mpdfc->getEvaluatorTypeUsed());
        }


        void test_DBPDF_delta12_click()
        {
            mpdfc->setEvaluatorType(OPTIMIZED);
            mpdfc->eval(mstru10);
            mpdfc->setDoubleAttr("delta1", mpdfc->getDoubleAttr("delta1"));
            mpdfc->setDoubleAttr("delta2", mpdfc->getDoubleAttr("delta2"));
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfc->getEvaluatorTypeUsed());
            mpdfc->setDoubleAttr("delta1", 0.5 + mpdfc->getDoubleAttr("delta1"));
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(BASIC, mpdfc->getEvaluatorTypeUsed());
            mpdfc->eval(mstru10);
            mpdfc->setDoubleAttr("delta2", 0.5 + mpdfc->getDoubleAttr("delta2"));
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(BASIC, mpdfc->getEvaluatorTypeUsed());
        }


        void test_DBPDF_SFTB_click()
        {
            mpdfc->setEvaluatorType(OPTIMIZED);
            mpdfc->eval(mstru10);
            mpdfc->setScatteringFactorTable(mpdfc->getScatteringFactorTable());
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfc->getEvaluatorTypeUsed());
            mpdfc->setScatteringFactorTableByType("neutron");
            mpdfc->eval(mstru10);
            TS_ASSERT_EQUALS(BASIC, mpdfc->getEvaluatorTypeUsed());
        }


};  // class TestDebyePDFCalculator

// End of file
