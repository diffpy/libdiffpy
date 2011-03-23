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
* class TestOverlapCalculator -- unit tests for the OverlapCalculator class
*
* $Id$
*
*****************************************************************************/

#include <memory>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/OverlapCalculator.hpp>
#include <diffpy/srreal/VR3Structure.hpp>
#include <diffpy/PythonInterface.hpp>
#include <diffpy/srreal/PythonStructureAdapter.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

class TestOverlapCalculator : public CxxTest::TestSuite
{
    private:

        boost::shared_ptr<OverlapCalculator> molc;
        VR3Structure memptystru;
        double meps;

    public:

        void setUp()
        {
            meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
            molc.reset(new OverlapCalculator);
            molc->getAtomRadiiTable()->setCustom("", 1.0);
        }


        void test_lineNoTouch()
        {
            VR3Structure stru;
            stru.push_back(R3::Vector(0.0, 0.0, 0.0));
            stru.push_back(R3::Vector(0.0, 0.0, 2.0));
            molc->eval(stru);
            TS_ASSERT_EQUALS(0.0, molc->totalSquareOverlap());
            QuantityType sqolps = molc->siteSquareOverlaps();
            TS_ASSERT_EQUALS(2u, sqolps.size());
            TS_ASSERT_EQUALS(0.0, sqolps[0]);
            TS_ASSERT_EQUALS(0.0, sqolps[1]);
            TS_ASSERT_EQUALS(0.0, molc->totalFlipDiff(0, 0));
            TS_ASSERT_EQUALS(0.0, molc->totalFlipDiff(0, 1));
            std::vector<R3::Vector> g = molc->gradients();
            TS_ASSERT_EQUALS(2u, g.size());
            TS_ASSERT_EQUALS(0.0, R3::norm(g[0]));
            TS_ASSERT_EQUALS(0.0, R3::norm(g[1]));
        }


        void test_lineTouch()
        {
            VR3Structure stru;
            stru.push_back(R3::Vector(0.0, 0.0, 0.0));
            stru.push_back(R3::Vector(0.0, 0.0, 1.5));
            molc->eval(stru);
            TS_ASSERT_EQUALS(0.25, molc->totalSquareOverlap());
            QuantityType sqolps = molc->siteSquareOverlaps();
            TS_ASSERT_EQUALS(2u, sqolps.size());
            TS_ASSERT_EQUALS(0.125, sqolps[0]);
            TS_ASSERT_EQUALS(0.125, sqolps[1]);
            TS_ASSERT_EQUALS(0.0, molc->totalFlipDiff(0, 1));
            std::vector<R3::Vector> g = molc->gradients();
            TS_ASSERT_EQUALS(2u, g.size());
            R3::Vector g0(0.0, 0.0, +1);
            TS_ASSERT_EQUALS(0.0, R3::distance(g0, g[0]));
            R3::Vector g1(0.0, 0.0, -1);
            TS_ASSERT_EQUALS(0.0, R3::distance(g1, g[1]));
        }


        void test_gradients()
        {
            VR3Structure stru;
            stru.push_back(R3::Vector(0.0, 0.0, 0.0));
            stru.push_back(R3::Vector(0.17, 0.31, 0.37));
            molc->eval(stru);
            std::vector<R3::Vector> g = molc->gradients();
            double dx = 1e-6;
            double tsqo0 = molc->totalSquareOverlap();
            R3::Vector r1 = stru[1];
            R3::Vector gn1;
            for (int k = 0; k != R3::Ndim; ++k)
            {
                stru[1] = r1;
                stru[1][k] += dx;
                molc->eval(stru);
                gn1[k] = (molc->totalSquareOverlap() - tsqo0) / dx;
            }
            TS_ASSERT(R3::norm(gn1) > 1.0);
            TS_ASSERT_DELTA(0.0, R3::norm(gn1 - g[1]), 1e-5);
            TS_ASSERT_DELTA(0.0, R3::norm(g[0] + g[1]), 1e-5);
        }


        void test_bccTouch()
        {
            using namespace boost;
            diffpy::initializePython();
            python::object Structure =
                diffpy::importFromPyModule("diffpy.Structure", "Structure");
            python::object bcc = Structure();
            bcc.attr("lattice").attr("setLatPar")(2.0, 2.0, 2.0);
            bcc.attr("addNewAtom")("C", python::make_tuple(0.0, 0.0, 0.0));
            bcc.attr("addNewAtom")("C", python::make_tuple(0.5, 0.5, 0.5));
            molc->getAtomRadiiTable()->setCustom("C", 1.0);
            molc->eval(bcc);
            TS_ASSERT_EQUALS(2, molc->getStructure()->countSites());
            TS_ASSERT_EQUALS(8 * pow(2.0 - sqrt(3.0), 2),
                    molc->totalSquareOverlap());
            // gradient is zero due to bcc symmetry
            std::vector<R3::Vector> g = molc->gradients();
            TS_ASSERT_DELTA(0.0, R3::norm(g[0]), meps);
            TS_ASSERT_DELTA(0.0, R3::norm(g[1]), meps);
            // move the second atom so it touches only one ball
            bcc[1].attr("xyz")[0] = -0.25;
            bcc[1].attr("xyz")[1] = 0.0;
            bcc[1].attr("xyz")[2] = 0.0;
            molc->eval(bcc);
            double tsqo = 1.5 * 1.5 + 0.5 * 0.5;
            TS_ASSERT_DELTA(tsqo, molc->totalSquareOverlap(), meps);
            g = molc->gradients();
            TS_ASSERT_EQUALS(2u, g.size());
            double gx = 1.5 + 0.5;
            TS_ASSERT_DELTA(-gx, g[0][0], meps);
            TS_ASSERT_DELTA(0.0, g[0][1], meps);
            TS_ASSERT_DELTA(0.0, g[0][2], meps);
            TS_ASSERT_DELTA(+gx, g[1][0], meps);
            TS_ASSERT_DELTA(0.0, g[1][1], meps);
            TS_ASSERT_DELTA(0.0, g[1][2], meps);
        };


        void test_serialization()
        {
            // build customized PDFCalculator
            molc->getAtomRadiiTable()->setCustom("", 1.0);
            VR3Structure stru;
            stru.push_back(R3::Vector(0.0, 0.0, 0.0));
            stru.push_back(R3::Vector(0.0, 0.0, 1.5));
            molc->eval(stru);
            boost::shared_ptr<OverlapCalculator> olc1;
            olc1 = dumpandload(molc);
            TS_ASSERT_DIFFERS(molc.get(), olc1.get());
            TS_ASSERT_EQUALS(2, olc1->getStructure()->countSites());
            TS_ASSERT_EQUALS(0.25, olc1->totalSquareOverlap());
            QuantityType sqolps = molc->siteSquareOverlaps();
            TS_ASSERT_EQUALS(2u, olc1->siteSquareOverlaps().size());
            TS_ASSERT_EQUALS(0.125, olc1->siteSquareOverlaps()[0]);
            TS_ASSERT_EQUALS(0.125, olc1->siteSquareOverlaps()[1]);
        }

};  // class TestOverlapCalculator

// End of file
