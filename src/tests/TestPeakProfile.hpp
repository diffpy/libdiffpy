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
* class TestPeakProfile -- unit tests for various PeakProfile classes
*
*****************************************************************************/

#include <cmath>
#include <stdexcept>
#include <cxxtest/TestSuite.h>

#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/PeakProfile.hpp>
#include "serialization_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;


class TestPeakProfile : public CxxTest::TestSuite
{

    private:

        PeakProfilePtr mpkgauss;
        double meps;

    public:

        void setUp()
        {
            mpkgauss = PeakProfile::createByType("gaussian");
            meps = 10 * diffpy::mathutils::DOUBLE_EPS;
        }


        void test_factory()
        {
            TS_ASSERT_THROWS(PeakProfile::createByType("invalid"),
                    invalid_argument);
            TS_ASSERT_EQUALS(string("gaussian"), mpkgauss->type());
        }


        void test_y()
        {
            double Afwhm1 = 2 * sqrt(M_LN2 / M_PI);
            const PeakProfile& pkgauss = *mpkgauss;
            TS_ASSERT_DELTA(Afwhm1, pkgauss(0, 1), meps);
            TS_ASSERT_DELTA(Afwhm1 / 2, pkgauss(1, 1), meps);
            TS_ASSERT_DELTA(Afwhm1 / 2, pkgauss(-1, 1), meps);
        }


        void test_xboundlo()
        {
            double epsy = 1e-8;
            mpkgauss->setPrecision(epsy);
            double xblo1 = mpkgauss->xboundlo(1);
            double xblo3 = mpkgauss->xboundlo(3);
            TS_ASSERT_DELTA(epsy, (*mpkgauss)(xblo1, 1), meps);
            TS_ASSERT_DELTA(epsy, (*mpkgauss)(xblo3, 3), meps);
            TS_ASSERT(xblo1 < 0);
            TS_ASSERT(xblo3 < 0);
            mpkgauss->setPrecision(10);
            TS_ASSERT_EQUALS(0.0, mpkgauss->xboundlo(1));
            mpkgauss->setPrecision(0.1);
            TS_ASSERT_EQUALS(0.0, mpkgauss->xboundlo(0));
            TS_ASSERT_EQUALS(0.0, mpkgauss->xboundlo(-1));
        }


        void test_xboundhi()
        {
            double epsy = 1e-8;
            mpkgauss->setPrecision(epsy);
            double xbhi1 = mpkgauss->xboundhi(1);
            double xbhi3 = mpkgauss->xboundhi(3);
            TS_ASSERT_DELTA(epsy, (*mpkgauss)(xbhi1, 1), meps);
            TS_ASSERT_DELTA(epsy, (*mpkgauss)(xbhi3, 3), meps);
            TS_ASSERT(xbhi1 > 0);
            TS_ASSERT(xbhi3 > 0);
        }


        void test_setPrecision()
        {
            using diffpy::mathutils::DOUBLE_EPS;
            double epsy = 1e-7;
            TS_ASSERT_DELTA(0.0, mpkgauss->getPrecision(), 10 * DOUBLE_EPS);
            mpkgauss->setPrecision(epsy);
            TS_ASSERT_EQUALS(epsy, mpkgauss->getPrecision());
            double xbhi1 = mpkgauss->xboundhi(1);
            mpkgauss->setPrecision(1e-4);
            TS_ASSERT(xbhi1 != mpkgauss->xboundhi(1));
        }


        void test_serialization()
        {
            mpkgauss->setPrecision(0.0123);
            PeakProfilePtr pk1 = dumpandload(mpkgauss);
            TS_ASSERT(pk1.get());
            TS_ASSERT_DIFFERS(mpkgauss.get(), pk1.get());
            TS_ASSERT_EQUALS(string("gaussian"), pk1->type());
            TS_ASSERT_EQUALS(0.0123, pk1->getPrecision());
        }

};  // class TestPeakProfile

// End of file
