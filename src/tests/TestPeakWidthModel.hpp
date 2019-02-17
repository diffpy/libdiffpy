/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2019 Brookhaven Science Associates,
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
* class TestPeakWidthModel -- test PeakWidthModel and derived classes
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/srreal/DebyeWallerPeakWidth.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

const double tuiso = 0.004;

//////////////////////////////////////////////////////////////////////////////
// class TestPeakWidthModel
//////////////////////////////////////////////////////////////////////////////

class TestPeakWidthModel : public CxxTest::TestSuite
{
    private:

        // data
        PeakWidthModelPtr mdwpw;
        PeakWidthModelPtr mjepw;
        JeongPeakWidth* mcjepw;
        double meps;

        // methods
        static AtomicStructureAdapterPtr dimer()
        {
            AtomicStructureAdapterPtr adpt(new AtomicStructureAdapter);
            Atom a;
            a.atomtype = "C";
            a.uij_cartn = R3::identity() * tuiso;
            adpt->append(a);
            a.xyz_cartn[2] = 1;
            adpt->append(a);
            return adpt;
        }

    public:

        void setUp()
        {
            mdwpw = PeakWidthModel::createByType("debye-waller");
            mjepw = PeakWidthModel::createByType("jeong");
            mcjepw = static_cast<JeongPeakWidth*>(mjepw.get());
            meps = 10 * diffpy::mathutils::DOUBLE_EPS;
        }


        void test_clone()
        {
            TS_ASSERT_EQUALS(mdwpw->type(), mdwpw->clone()->type());
            TS_ASSERT_EQUALS(mjepw->type(), mjepw->clone()->type());
        }


        void test_attributes()
        {
            TS_ASSERT_EQUALS(0, mdwpw->namesOfDoubleAttributes().size());
            TS_ASSERT_EQUALS(3, mjepw->namesOfDoubleAttributes().size());
            mcjepw->setDelta1(0.1);
            mcjepw->setDelta2(0.2);
            mcjepw->setQbroad(0.3);
            TS_ASSERT_EQUALS(0.1, mjepw->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(0.2, mjepw->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(0.3, mjepw->getDoubleAttr("qbroad"));
        }


        void test_calculate()
        {
            using diffpy::mathutils::GAUSS_SIGMA_TO_FWHM;
            AtomicStructureAdapterPtr adpt = dimer();
            BaseBondGeneratorPtr bnds = adpt->createBondGenerator();
            bnds->setRmax(2.0);
            bnds->rewind();
            const double w = GAUSS_SIGMA_TO_FWHM * sqrt(2 * tuiso);
            TS_ASSERT_DELTA(w, mdwpw->calculate(*bnds), meps);
            TS_ASSERT_DELTA(w, mjepw->calculate(*bnds), meps);
            mcjepw->setDelta2(0.75);
            TS_ASSERT_DELTA(0.5 * w, mjepw->calculate(*bnds), meps);
            mcjepw->setDelta2(0);
            mcjepw->setQbroad(sqrt(3));
            TS_ASSERT_DELTA(2 * w, mjepw->calculate(*bnds), meps);
        }


        void test_serialization()
        {
            PeakWidthModelPtr pwm1 = dumpandload(mdwpw);
            DebyeWallerPeakWidth* cpwm1 =
                dynamic_cast<DebyeWallerPeakWidth*>(pwm1.get());
            TS_ASSERT(cpwm1);
            mcjepw->setDelta1(1.1);
            mcjepw->setDoubleAttr("delta2", 2.2);
            mcjepw->setDoubleAttr("qbroad", 0.03);
            PeakWidthModelPtr pwm2 = dumpandload(mjepw);
            JeongPeakWidth* cpwm2 = dynamic_cast<JeongPeakWidth*>(pwm2.get());
            TS_ASSERT(cpwm2);
            TS_ASSERT_EQUALS(1.1, pwm2->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(2.2, pwm2->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(0.03, pwm2->getDoubleAttr("qbroad"));
        }

};  // class TestPeakWidthModel

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestPeakWidthModel;

// End of file
