/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestPQEvaluator -- unit tests for the PairQuantity evaluator class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <functional>
#include <boost/make_shared.hpp>

#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/PairCounter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/OverlapCalculator.hpp>
#include "test_helpers.hpp"

namespace diffpy {
namespace srreal {

using namespace std;
using diffpy::mathutils::EpsilonEqual;

// calculator with invalid support for OPTIMIZED evaluation

class BadPairCounter : public PairCounter
{
    protected:

        // intentionally faulty support for PQEvaluatorOptimized
        virtual void stashPartialValue()  { }
        virtual void restorePartialValue()  { }
};

//////////////////////////////////////////////////////////////////////////////
// class TestPQEvaluator
//////////////////////////////////////////////////////////////////////////////

class TestPQEvaluator : public CxxTest::TestSuite
{
    private:

        // data
        EpsilonEqual allclose;
        PDFCalculator mpdfcb;
        PDFCalculator mpdfco;
        QuantityType mzeros;
        AtomicStructureAdapterPtr mstru10;
        AtomicStructureAdapterPtr mstru10d1;
        AtomicStructureAdapterPtr mstru10r;
        AtomicStructureAdapterPtr mstru9;

        // methods
        QuantityType pdfcdiff(StructureAdapterPtr stru)
        {
            mpdfcb.eval(stru);
            mpdfco.eval(stru);
            QuantityType gb = mpdfcb.getPDF();
            QuantityType go = mpdfco.getPDF();
            assert(mpdfcb.getRgrid() == mpdfco.getRgrid());
            QuantityType rv(gb.size());
            transform(gb.begin(), gb.end(), go.begin(),
                    rv.begin(), minus<double>());
            return rv;
        }

     public:

        void setUp()
        {
            // setup PDFCalculator instances and the zero vector for comparison
            mpdfcb.setEvaluatorType(BASIC);
            mpdfco.setEvaluatorType(OPTIMIZED);
            mzeros.assign(mpdfcb.getRgrid().size(), 0.0);
            // setup structure instances
            const int SZ = 10;
            mstru10 = boost::make_shared<AtomicStructureAdapter>();
            Atom ai;
            ai.atomtype = "C";
            ai.uij_cartn = R3::identity();
            ai.uij_cartn(0, 0) = ai.uij_cartn(1, 1) =
                ai.uij_cartn(2, 2) = 0.004;
            for (int i = 0; i < SZ; ++i)
            {
                ai.xyz_cartn[0] = i;
                mstru10->append(ai);
            }
            mstru10d1 = boost::make_shared<AtomicStructureAdapter>(*mstru10);
            (*mstru10d1)[0].atomtype = "Au";
            mstru10r = boost::make_shared<AtomicStructureAdapter>();
            mstru10r->assign(mstru10->rbegin(), mstru10->rend());
            mstru9 = boost::make_shared<AtomicStructureAdapter>(*mstru10);
            mstru9->erase(9);
        }


        void test_initial_evaluator_type_used()
        {
            PDFCalculator pdfc;
            TS_ASSERT_EQUALS(NONE, pdfc.getEvaluatorTypeUsed());
            pdfc.setEvaluatorType(OPTIMIZED);
            TS_ASSERT_EQUALS(NONE, pdfc.getEvaluatorTypeUsed());
            pdfc.setEvaluatorType(CHECK);
            TS_ASSERT_EQUALS(NONE, pdfc.getEvaluatorTypeUsed());
            pdfc.setEvaluatorType(BASIC);
            TS_ASSERT_EQUALS(NONE, pdfc.getEvaluatorTypeUsed());
        }


        void test_PDF_change_atom()
        {
            using std::placeholders::_1;
            TS_ASSERT_EQUALS(BASIC, mpdfcb.getEvaluatorType());
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorType());
            TS_ASSERT_EQUALS(mzeros, this->pdfcdiff(mstru10));
            // first call of mpdfco should use the BASIC evaluation
            TS_ASSERT_EQUALS(BASIC, mpdfcb.getEvaluatorTypeUsed());
            TS_ASSERT_EQUALS(BASIC, mpdfco.getEvaluatorTypeUsed());
            // verify there are some non-zero values in the PDF.
            QuantityType gb = mpdfcb.getPDF();
            TS_ASSERT(!gb.empty());
            int cnonzero = count_if(gb.begin(), gb.end(),
                    bind(not_equal_to<double>(), _1, 0.0));
            TS_ASSERT(cnonzero);
            // test second call on the same structure
            mpdfco.eval(mstru10);
            TS_ASSERT_EQUALS(gb, mpdfco.getPDF());
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            // test structure with one different atom
            TS_ASSERT(allclose(mzeros, this->pdfcdiff(mstru10d1)));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb1 = mpdfcb.getPDF();
            TS_ASSERT(!allclose(gb, gb1));
            // change position of 1 atom
            mstru10d1->at(0).xyz_cartn[1] = 0.5;
            TS_ASSERT(allclose(mzeros, this->pdfcdiff(mstru10d1)));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb2 = mpdfcb.getPDF();
            TS_ASSERT(!allclose(gb1, gb2));
        }


        void test_PDF_reverse_atoms()
        {
            mpdfco.eval(mstru10);
            TS_ASSERT(allclose(mzeros, this->pdfcdiff(mstru10r)));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
        }


        void test_PDF_remove_atom()
        {
            mpdfco.eval(mstru10);
            TS_ASSERT(allclose(mzeros, this->pdfcdiff(mstru9)));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
        }


        void test_PDF_type_mask()
        {
            mpdfcb.setTypeMask("O2-", "all", false);
            mpdfco.setTypeMask("O2-", "all", false);
            PeriodicStructureAdapterPtr litao =
                boost::dynamic_pointer_cast<PeriodicStructureAdapter>(
                        loadTestPeriodicStructure("LiTaO3.stru"));
            TS_ASSERT_EQUALS(mzeros, this->pdfcdiff(litao));
            TS_ASSERT_EQUALS(BASIC, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb0 = mpdfcb.getPDF();
            // change masked-away oxygen
            Atom& o29 = litao->at(29);
            TS_ASSERT_EQUALS("O2-", o29.atomtype);
            o29.xyz_cartn = R3::Vector(0.1, 0.2, 0.3);
            TS_ASSERT_EQUALS(mzeros, this->pdfcdiff(litao));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb1 = mpdfcb.getPDF();
            TS_ASSERT_EQUALS(gb0, gb1);
            // change active lithium
            Atom& li0 = litao->at(0);
            li0.occupancy = 0.1;
            TS_ASSERT(allclose(mzeros, this->pdfcdiff(litao)));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb2 = mpdfcb.getPDF();
            TS_ASSERT(!allclose(gb0, gb2));
        }


        void test_PDF_pair_mask()
        {
            mpdfcb.setPairMask(0, 3, false);
            mpdfco.setPairMask(0, 3, false);
            TS_ASSERT_EQUALS(mzeros, this->pdfcdiff(mstru10));
            TS_ASSERT_EQUALS(BASIC, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb0 = mpdfcb.getPDF();
            // check structure with one changed atom
            TS_ASSERT(allclose(mzeros, this->pdfcdiff(mstru10d1)));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb1 = mpdfcb.getPDF();
            TS_ASSERT(!allclose(gb0, gb1));
            // check structure with one removed atom
            TS_ASSERT(allclose(mzeros, this->pdfcdiff(mstru9)));
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb2 = mpdfcb.getPDF();
            TS_ASSERT(!allclose(gb1, gb2));
            // check full evaluation for structure with reversed atoms
            mpdfco.eval(mstru10);
            TS_ASSERT_EQUALS(mzeros, this->pdfcdiff(mstru10r));
            TS_ASSERT_EQUALS(BASIC, mpdfco.getEvaluatorTypeUsed());
            QuantityType gb3 = mpdfcb.getPDF();
            TS_ASSERT(allclose(gb0, gb3));
            // check effect of pair mask updates
            mpdfco.eval(mstru10r);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            mpdfco.setPairMask(0, 3, false);
            mpdfco.eval(mstru10r);
            TS_ASSERT_EQUALS(OPTIMIZED, mpdfco.getEvaluatorTypeUsed());
            mpdfco.setPairMask(0, 3, true);
            mpdfco.eval(mstru10r);
            TS_ASSERT_EQUALS(BASIC, mpdfco.getEvaluatorTypeUsed());
        }


        void test_optimized_supported()
        {
            mpdfcb.eval(mstru10);
            TS_ASSERT_EQUALS(BASIC, mpdfcb.getEvaluatorTypeUsed());
            mpdfcb.setEvaluatorType(OPTIMIZED);
            TS_ASSERT_EQUALS(BASIC, mpdfcb.getEvaluatorTypeUsed());
        }


        void test_optimized_unsupported()
        {
            OverlapCalculator olc;
            TS_ASSERT_EQUALS(BASIC, olc.getEvaluatorType());
            TS_ASSERT_THROWS(
                    olc.setEvaluatorType(OPTIMIZED), invalid_argument);
            TS_ASSERT_EQUALS(BASIC, olc.getEvaluatorType());
        }


        void test_checkevaluator_unsupported()
        {
            OverlapCalculator olc;
            TS_ASSERT_EQUALS(BASIC, olc.getEvaluatorType());
            TS_ASSERT_THROWS(
                    olc.setEvaluatorType(CHECK), invalid_argument);
            TS_ASSERT_EQUALS(BASIC, olc.getEvaluatorType());
        }


        void test_checkevaluator_passes()
        {
            mpdfco.setEvaluatorType(CHECK);
            TS_ASSERT_EQUALS(CHECK, mpdfco.getEvaluatorType());
            mpdfco.eval(mstru10);
            TS_ASSERT_EQUALS(BASIC, mpdfco.getEvaluatorTypeUsed());
            mpdfco.eval(mstru10);
            TS_ASSERT_EQUALS(CHECK, mpdfco.getEvaluatorTypeUsed());
            mpdfco.eval(mstru9);
            TS_ASSERT_EQUALS(CHECK, mpdfco.getEvaluatorTypeUsed());
            mpdfco.eval(emptyStructureAdapter());
            TS_ASSERT_EQUALS(BASIC, mpdfco.getEvaluatorTypeUsed());
        }


        void test_checkevaluator_fails()
        {
            BadPairCounter badcounter;
            badcounter.setEvaluatorType(CHECK);
            TS_ASSERT_EQUALS(45, badcounter(mstru10));
            TS_ASSERT_EQUALS(BASIC, badcounter.getEvaluatorTypeUsed());
            TS_ASSERT_THROWS(badcounter(mstru10), logic_error);
            TS_ASSERT_EQUALS(CHECK, badcounter.getEvaluatorTypeUsed());
        }

};  // class TestPQEvaluator

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestPQEvaluator;

// End of file
