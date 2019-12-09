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
#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/srreal/DebyeWallerPeakWidth.hpp>
#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include "serialization_helpers.hpp"

namespace diffpy {
namespace srreal {

const double tuiso = 0.004;
const double UtoB = 8 * M_PI * M_PI;

//////////////////////////////////////////////////////////////////////////////
// class TestPeakWidthModel
//////////////////////////////////////////////////////////////////////////////

class TestPeakWidthModel : public CxxTest::TestSuite
{
    private:

        // data
        PeakWidthModelPtr mcnpw;
        ConstantPeakWidth* mccnpw;
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
            mcnpw = PeakWidthModel::createByType("constant");
            mccnpw = static_cast<ConstantPeakWidth*>(mcnpw.get());
            mdwpw = PeakWidthModel::createByType("debye-waller");
            mjepw = PeakWidthModel::createByType("jeong");
            mcjepw = static_cast<JeongPeakWidth*>(mjepw.get());
            meps = 10 * diffpy::mathutils::DOUBLE_EPS;
        }


        void test_clone()
        {
            TS_ASSERT_EQUALS("constant", mcnpw->clone()->type());
            TS_ASSERT_EQUALS("debye-waller", mdwpw->clone()->type());
            TS_ASSERT_EQUALS("jeong", mjepw->clone()->type());
        }


        void test_attributes()
        {
            TS_ASSERT_EQUALS(0, mdwpw->namesOfDoubleAttributes().size());
            TS_ASSERT_EQUALS(4, mjepw->namesOfDoubleAttributes().size());
            mcjepw->setDelta1(0.1);
            mcjepw->setDelta2(0.2);
            mcjepw->setQbroad(0.3);
            mcjepw->setQbroad_seperable(0.4);
            TS_ASSERT_EQUALS(0.1, mjepw->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(0.2, mjepw->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(0.3, mjepw->getDoubleAttr("qbroad"));
            TS_ASSERT_EQUALS(0.4, mjepw->getDoubleAttr("qbroad_seperable"));
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
            mcjepw->setQbroad_seperable(0);
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
            mcjepw->setDoubleAttr("qbroad_seperable", 0.003);
            PeakWidthModelPtr pwm2 = dumpandload(mjepw);
            JeongPeakWidth* cpwm2 = dynamic_cast<JeongPeakWidth*>(pwm2.get());
            TS_ASSERT(cpwm2);
            TS_ASSERT_EQUALS(1.1, pwm2->getDoubleAttr("delta1"));
            TS_ASSERT_EQUALS(2.2, pwm2->getDoubleAttr("delta2"));
            TS_ASSERT_EQUALS(0.03, pwm2->getDoubleAttr("qbroad"));
            TS_ASSERT_EQUALS(0.003, pwm2->getDoubleAttr("qbroad_seperable"));
        }


        void test_CPWM_attributes()
        {
            TS_ASSERT_EQUALS("constant", mcnpw->type());
            TS_ASSERT(mcnpw->hasDoubleAttr("width"));
            TS_ASSERT(mcnpw->hasDoubleAttr("bisowidth"));
            TS_ASSERT(mcnpw->hasDoubleAttr("uisowidth"));
            TS_ASSERT_EQUALS(3, mcnpw->namesOfDoubleAttributes().size());
            mccnpw->setWidth(0.1);
            TS_ASSERT_EQUALS(0.1, mcnpw->getDoubleAttr("width"));
            TS_ASSERT_EQUALS(0.1, mcnpw->maxWidth(dimer(), 0, 4));
            const double bw0 = mcnpw->getDoubleAttr("bisowidth");
            const double uw0 = mcnpw->getDoubleAttr("uisowidth");
            mccnpw->setWidth(2 * mccnpw->getWidth());
            TS_ASSERT_DELTA(bw0, UtoB * uw0, meps);
            TS_ASSERT_DELTA(4 * bw0, mcnpw->getDoubleAttr("bisowidth"), meps);
            TS_ASSERT_DELTA(4 * uw0, mcnpw->getDoubleAttr("uisowidth"), meps);
        }


        void test_CPWM_bisowidth()
        {
            const double biso = UtoB * tuiso;
            mcnpw->setDoubleAttr("bisowidth", biso);
            TS_ASSERT_DELTA(tuiso, mcnpw->getDoubleAttr("uisowidth"), meps);
            mcnpw->setDoubleAttr("uisowidth", -tuiso);
            TS_ASSERT_DELTA(-biso, mcnpw->getDoubleAttr("bisowidth"), meps);
        }


        void test_CPWM_uisowidth()
        {
            StructureAdapterPtr adpt = dimer();
            BaseBondGeneratorPtr bnds = adpt->createBondGenerator();
            PeakWidthModelPtr dwpw(new DebyeWallerPeakWidth);
            bnds->setRmax(2.0);
            bnds->rewind();
            const double w0 = dwpw->calculate(*bnds);
            TS_ASSERT_EQUALS(0.0, mcnpw->calculate(*bnds));
            mcnpw->setDoubleAttr("uisowidth", tuiso);
            TS_ASSERT_DELTA(w0, mcnpw->calculate(*bnds), meps);
            mcnpw->setDoubleAttr("uisowidth", -tuiso);
            TS_ASSERT_DELTA(-tuiso, mcnpw->getDoubleAttr("uisowidth"), meps);
            TS_ASSERT_DELTA(-w0, mcnpw->getDoubleAttr("width"), meps);
        }


        void test_CPWM_serialization()
        {
            PeakWidthModelPtr pwm1 = dumpandload(mcnpw);
            TS_ASSERT_EQUALS("constant", pwm1->type());
            TS_ASSERT_EQUALS(0.0, pwm1->getDoubleAttr("width"));
            mcnpw->setDoubleAttr("width", 0.1);
            double uw0 = mcnpw->getDoubleAttr("uisowidth");
            PeakWidthModelPtr pwm2 = dumpandload(mcnpw);
            TS_ASSERT_EQUALS(0.1, pwm2->getDoubleAttr("width"));
            TS_ASSERT_DELTA(uw0, pwm2->getDoubleAttr("uisowidth"), meps);
        }



};  // class TestPeakWidthModel

}   // namespace srreal
}   // namespace diffpy

using diffpy::srreal::TestPeakWidthModel;

// End of file
