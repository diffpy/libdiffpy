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
* class TestBVSCalculator -- unit tests suite
*
* $Id$
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/PythonInterface.hpp>
#include <diffpy/srreal/PythonStructureAdapter.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>
#include "tests_dir.hpp"

using namespace std;
using namespace boost;
using diffpy::srreal::BVSCalculator;

// Local Helpers -------------------------------------------------------------

namespace {

python::object loadTestStructure(const string& tailname)
{
    python::object Structure_class =
        diffpy::importFromPyModule("diffpy.Structure", "Structure");
    python::object stru = Structure_class();
    string fp = prepend_testdata_dir(tailname);
    stru.attr("read")(fp);
    return stru;
}

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// class TestBVSCalculator
//////////////////////////////////////////////////////////////////////////////

class TestBVSCalculator : public CxxTest::TestSuite
{
private:

    python::object mnacl;
    auto_ptr<BVSCalculator> mbvc;

public:

    void setUp()
    {
        diffpy::initializePython();
        if (mnacl.ptr() == Py_None)
        {
            mnacl = loadTestStructure("NaCl.cif");
        }
        mbvc.reset(new BVSCalculator);
    }


    void test_NaCl()
    {
        const double eps = 0.02;
        mbvc->eval(mnacl);
        TS_ASSERT_EQUALS(8u, mbvc->value().size());
        TS_ASSERT_DELTA(+1, mbvc->value()[0], eps);
        TS_ASSERT_DELTA(+1, mbvc->value()[1], eps);
        TS_ASSERT_DELTA(+1, mbvc->value()[2], eps);
        TS_ASSERT_DELTA(+1, mbvc->value()[3], eps);
        TS_ASSERT_DELTA(-1, mbvc->value()[4], eps);
        TS_ASSERT_DELTA(-1, mbvc->value()[5], eps);
        TS_ASSERT_DELTA(-1, mbvc->value()[6], eps);
        TS_ASSERT_DELTA(-1, mbvc->value()[7], eps);
    }

};  // class TestBVSCalculator

// End of file
