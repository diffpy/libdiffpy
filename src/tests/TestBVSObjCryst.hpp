/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestBVSObjCryst -- unit tests for BVS calculation for an ObjCryst
*   crystal structure
*
* $Id$
*
*****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>
#include <ObjCryst/ObjCryst/CIF.h>

#include "tests_dir.hpp"

using namespace std;
using namespace diffpy::srreal;

// Local Helpers -------------------------------------------------------------

namespace {

ObjCryst::Crystal* loadTestCIF(const string& tailname)
{
    string fp = prepend_testdata_dir(tailname);
    ifstream in(fp.c_str());
    ObjCryst::CIF cif(in);
    // redirect cout to hide the CreateCrystalFromCIF chatty output
    streambuf* cout_sb = cout.rdbuf();
    ostringstream cifout;
    cout.rdbuf(cifout.rdbuf());
    ObjCryst::Crystal* crst = ObjCryst::CreateCrystalFromCIF(cif, false);
    // restore cout
    cout.rdbuf(cout_sb);
    return crst;
}

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// class Test
//////////////////////////////////////////////////////////////////////////////

class TestBVSObjCryst : public CxxTest::TestSuite
{
    private:

        auto_ptr<ObjCryst::Crystal> mcrystal;
        auto_ptr<BVSCalculator> mbvc;

    public:

        void setUp()
        {
            if (!mcrystal.get())
            {
                mcrystal.reset(loadTestCIF("NaCl.cif"));
            }
            mbvc.reset(new BVSCalculator);
        }


        void test_NaCl()
        {
            const double eps = 0.02;
            mbvc->eval(*mcrystal);
            TS_ASSERT_EQUALS(2u, mbvc->value().size());
            TS_ASSERT_DELTA(+1, mbvc->value()[0], eps);
            TS_ASSERT_DELTA(-1, mbvc->value()[1], eps);
        }

};  // class TestBVSObjCryst

// End of file
