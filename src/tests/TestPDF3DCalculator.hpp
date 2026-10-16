/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class TestPDF3DCalculator -- unit tests for PDF3DCalculator
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/PDF3DCalculator.hpp>

using namespace diffpy::srreal;


class TestPDF3DCalculator : public CxxTest::TestSuite
{

public:

    void test_grid_step_stays_synchronized()
    {
        PDF3DCalculator calc;
        TS_ASSERT_EQUALS(calc.getGridStep(), calc.getRstep());

        calc.setGridStep(0.25);
        TS_ASSERT_EQUALS(0.25, calc.getGridStep());
        TS_ASSERT_EQUALS(calc.getGridStep(), calc.getRstep());

        PDFCalculator& base = calc;
        base.setRstep(0.2);
        TS_ASSERT_EQUALS(0.2, calc.getGridStep());
        TS_ASSERT_EQUALS(calc.getGridStep(), calc.getRstep());
        TS_ASSERT_THROWS(calc.setGridStep(0.0), std::invalid_argument);
    }


    void test_modes_validation()
    {
        PDF3DCalculator calc;
        calc.setCalculationMode3D(1);
        calc.setHistogramWeightMode3D(1);
        TS_ASSERT_EQUALS(1, calc.getCalculationMode3D());
        TS_ASSERT_EQUALS(1, calc.getHistogramWeightMode3D());
        TS_ASSERT_THROWS(calc.setCalculationMode3D(2), std::invalid_argument);
        TS_ASSERT_THROWS(calc.setHistogramWeightMode3D(-1),
                std::invalid_argument);
    }

};  // class TestPDF3DCalculator

// End of file
