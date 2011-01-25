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

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>

#include "objcryst_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

//////////////////////////////////////////////////////////////////////////////
// class Test
//////////////////////////////////////////////////////////////////////////////

class TestBVSObjCryst : public CxxTest::TestSuite
{
    private:

        auto_ptr<ObjCryst::Crystal> mnacl;
        auto_ptr<BVSCalculator> mbvc;

    public:

        void setUp()
        {
            if (!mnacl.get())
            {
                mnacl.reset(loadTestCrystal("NaCl.cif"));
            }
            mbvc.reset(new BVSCalculator);
        }


        void test_NaCl()
        {
            const double eps = 1e-4;
            mbvc->eval(*mnacl);
            TS_ASSERT_EQUALS(2u, mbvc->value().size());
            TS_ASSERT_DELTA(+1.01352, mbvc->value()[0], eps);
            TS_ASSERT_DELTA(-1.01352, mbvc->value()[1], eps);
        }


        void test_NaCl_mixed()
        {
            auto_ptr<ObjCryst::Crystal> nacl_mixed;
            nacl_mixed.reset(loadTestCrystal("NaCl_mixed.cif"));
            mbvc->eval(*mnacl);
            BVSCalculator bvc;
            bvc.eval(*nacl_mixed);
            TS_ASSERT_DELTA(mbvc->bvrmsdiff(), bvc.bvrmsdiff(), 1e-12);
        }

};  // class TestBVSObjCryst

// End of file
